"""界面残基的选取，以及 Pose 的拆分与 repack —— 各档共用。

抽出来是因为 01/02/03 各自复制过一份，改一处要改三处；插入编号那个 bug 就是
重复维护的代价。

# 两套判据，回答两个不同的问题

    位点判据    哪些残基真的在界面上         接触问题 —— 重原子 4 Å
    repack 判据 算能量时放开哪些侧链去松弛    堆积问题 —— 邻居原子(CB) 8 Å

两个数值**不可换算**。同一结构的抗原侧实测：cb 4 Å 只选出 1 个残基（两个 CB
相距 4 Å 意味着主链快重叠了），cb 8 Å 是 15 个，重原子 4 Å 是 23 个，cb 10 Å
是 23 个。换判据必须同时换距离。

⚠️ Rosetta 自带的 CloseContactResidueSelector **把氢也算进去**，实测同一结构
4 Å 下它给 32 个、重原子口径只有 23 个，差 39%。所以重原子判据在这里自己实现。
"""

CB, HEAVY = 'cb', 'heavy'


def chains_selector(chains):
    """'HL' -> ChainSelector('H') OR ChainSelector('L')"""
    from pyrosetta.rosetta.core.select.residue_selector import (
        ChainSelector, OrResidueSelector)

    sel = ChainSelector(chains[0])
    for c in chains[1:]:
        sel = OrResidueSelector(sel, ChainSelector(c))
    return sel


def chain_indices(pose, chains):
    """某几条链的残基编号"""
    info = pose.pdb_info()
    return [i for i in range(1, pose.total_residue() + 1)
            if info.chain(i) in chains]


def heavy_contacts(pose, mine, other, cut):
    """mine 里与 other 存在重原子接触的残基

    Rosetta 的残基把重原子排在前面（1..nheavyatoms()），氢在后面。
    """
    hit = []
    for i in mine:
        ri = pose.residue(i)
        xi = ri.nbr_atom_xyz()
        for j in other:
            rj = pose.residue(j)
            # 任一原子都在自己 nbr_atom 的 nbr_radius 之内，据此安全地跳过远的残基对
            if xi.distance(rj.nbr_atom_xyz()) > cut + ri.nbr_radius() + rj.nbr_radius():
                continue
            if any(ri.xyz(a).distance(rj.xyz(b)) <= cut
                   for a in range(1, ri.nheavyatoms() + 1)
                   for b in range(1, rj.nheavyatoms() + 1)):
                hit.append(i)
                break
    return hit


def near_residues(pose, mine_chains, other_chains, dist, criterion):
    """mine 那一侧、离 other 足够近的残基编号"""
    from pyrosetta.rosetta.core.select import get_residues_from_subset
    from pyrosetta.rosetta.core.select.residue_selector import (
        AndResidueSelector, NeighborhoodResidueSelector)

    if criterion == HEAVY:
        return heavy_contacts(pose, chain_indices(pose, mine_chains),
                              chain_indices(pose, other_chains), dist)
    sel = AndResidueSelector(
        NeighborhoodResidueSelector(chains_selector(other_chains), dist, False),
        chains_selector(mine_chains))
    return sorted(get_residues_from_subset(sel.apply(pose)))


def interface_both(pose, spec, dist, criterion):
    """界面两侧的残基编号"""
    g1, g2 = spec.split('_')
    return sorted(set(near_residues(pose, g1, g2, dist, criterion)) |
                  set(near_residues(pose, g2, g1, dist, criterion)))


def one_side(pose, spec, dist, criterion, side):
    """只取一侧的界面残基。side 为 antibody 时取下划线左边"""
    g1, g2 = spec.split('_')
    mine, other = (g1, g2) if side == 'antibody' else (g2, g1)
    return near_residues(pose, mine, other, dist, criterion)


def group_pose(pose, chains):
    """把指定链号的残基拼成一个独立的 Pose，同时返回原编号映射"""
    from pyrosetta.rosetta.core.pose import append_pose_to_pose

    parts = pose.split_by_chain()
    info = pose.pdb_info()
    out, orig = None, []
    for k in range(1, pose.num_chains() + 1):
        if info.chain(pose.chain_begin(k)) not in chains:
            continue
        orig += list(range(pose.chain_begin(k), pose.chain_end(k) + 1))
        if out is None:
            out = parts[k].clone()
        else:
            append_pose_to_pose(out, parts[k], True)
    out.conformation().detect_disulfides()    # 拆开后二硫键记录失效，必须重建
    return out, orig


def repack(pose, scorefxn, subset):
    """只 repack 指定残基的侧链，其余一律固定"""
    from pyrosetta.rosetta.core.pack.task import TaskFactory
    from pyrosetta.rosetta.core.pack.task.operation import (
        IncludeCurrent, OperateOnResidueSubset, PreventRepackingRLT, RestrictToRepacking)
    from pyrosetta.rosetta.core.select.residue_selector import (
        NotResidueSelector, ResidueIndexSelector)
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover

    tf = TaskFactory()
    tf.push_back(RestrictToRepacking())                 # 不改序列
    tf.push_back(IncludeCurrent())                      # 当前构象也纳入候选，否则丢掉已有的最小化成果
    keep = ResidueIndexSelector(','.join(map(str, subset)))
    tf.push_back(OperateOnResidueSubset(PreventRepackingRLT(), NotResidueSelector(keep)))
    pack = PackRotamersMover(scorefxn)
    pack.task_factory(tf)
    pack.apply(pose)


def dG_from_pose(pose, spec, idx, scorefxn):
    """ΔG = E_complex − E_antibody − E_antigen，idx 是 repack 范围"""
    repack(pose, scorefxn, idx)
    e_complex = scorefxn(pose)

    g1, g2 = spec.split('_')
    iface = set(idx)
    pose_ab, orig_ab = group_pose(pose, g1)
    pose_ag, orig_ag = group_pose(pose, g2)
    repack(pose_ab, scorefxn,
           [j for j, o in enumerate(orig_ab, 1) if o in iface])
    repack(pose_ag, scorefxn,
           [j for j, o in enumerate(orig_ag, 1) if o in iface])
    return e_complex - scorefxn(pose_ab) - scorefxn(pose_ag)
