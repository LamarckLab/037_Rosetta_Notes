"""批量计算复合物的结合能 ΔG —— 不 relax，只 repack（Rosetta 官方推荐的标准做法）。

选界面残基 -> 结合态 repack -> 打分 -> 拆开两侧各自 repack -> 打分 -> 相减
官方文档对非 Rosetta 来源的结构建议 pack_input + pack_separated，但不要求先 relax
"""

# ==================== 配置（日常改这里） ====================

INPUTS    = '/data/lmk/rosetta_inputs'                          # 待评估的复合物 pdb 目录
OUTPUT    = '/data/lmk/rosetta_outputs/dG_repack_results.csv'    # 指标输出
INTERFACE = 'HL_A'      # 链分组，要与 pdb 实际链号一致
                        # ⚠️ 左边必须是抗体、右边是抗原，否则 CSV 里两侧的列名会对调
RADIUS    = 8.0         # 界面判定半径 A（CB 之间的距离），决定 repack 范围

# ===========================================================

import argparse
import csv
import glob
import os
import time

import pyrosetta
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import (
    IncludeCurrent, OperateOnResidueSubset, PreventRepackingRLT, RestrictToRepacking)
from pyrosetta.rosetta.core.pose import append_pose_to_pose
from pyrosetta.rosetta.core.select import get_residues_from_subset
from pyrosetta.rosetta.core.select.residue_selector import (
    AndResidueSelector, ChainSelector, NeighborhoodResidueSelector, NotResidueSelector,
    OrResidueSelector, ResidueIndexSelector)
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.protocols.docking import setup_foldtree
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.utility import vector1_int

FIELDS = ['pdb_id', 'nres', 'nres_repack', 'nres_antibody', 'nres_antigen',
          'E_complex(REU)', 'E_antibody(REU)', 'E_antigen(REU)',
          'dG(REU)', 'dG_IA(REU)', 'dSASA_int(A^2)',
          'repack_time(s)', 'analyze_time(s)']


def chains_selector(chains):
    """'HL' -> ChainSelector('H') OR ChainSelector('L')"""
    sel = ChainSelector(chains[0])
    for c in chains[1:]:
        sel = OrResidueSelector(sel, ChainSelector(c))
    return sel


def interface_residues(pose, spec, radius):
    """界面两侧的残基编号（Rosetta 编号）"""
    g1, g2 = spec.split('_')
    s1, s2 = chains_selector(g1), chains_selector(g2)
    side1 = AndResidueSelector(NeighborhoodResidueSelector(s2, radius, False), s1)
    side2 = AndResidueSelector(NeighborhoodResidueSelector(s1, radius, False), s2)
    return sorted(set(get_residues_from_subset(side1.apply(pose))) |
                  set(get_residues_from_subset(side2.apply(pose))))


def group_pose(pose, chains):
    """把指定链号的残基拼成一个独立的 Pose，同时返回原编号映射"""
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
    tf = TaskFactory()
    tf.push_back(RestrictToRepacking())                 # 不改序列
    tf.push_back(IncludeCurrent())                      # 把当前构象也纳入候选，否则丢掉 relax 的最小化成果
    keep = ResidueIndexSelector(','.join(map(str, subset)))
    tf.push_back(OperateOnResidueSubset(PreventRepackingRLT(), NotResidueSelector(keep)))    # 只放开 subset
    pack = PackRotamersMover(scorefxn)
    pack.task_factory(tf)
    pack.apply(pose)


def evaluate(path, spec, radius, scorefxn):
    """算一个结构，返回一行指标。并行时每个进程各自 init 并传入自己的 scorefxn"""
    pose = pyrosetta.pose_from_pdb(path)

    idx = interface_residues(pose, spec, radius)
    if not idx:
        raise ValueError(f'界面为空，检查链分组 {spec} 是否与该文件匹配')

    pose_ia = pose.clone()           # 留一份原样，供 IA 独立计算

    t0 = time.time()
    repack(pose, scorefxn, idx)      # 结合态 repack，与分离态对等
    e_complex = scorefxn(pose)
    t_repack = time.time() - t0

    g1, g2 = spec.split('_')
    iface = set(idx)
    pose_ab, orig_ab = group_pose(pose, g1)                          # 拆开
    pose_ag, orig_ag = group_pose(pose, g2)

    # 界面残基在子 Pose 里的新编号；只 repack 它们，与复合物侧 relax 的范围保持一致
    sub_ab = [j for j, o in enumerate(orig_ab, 1) if o in iface]
    sub_ag = [j for j, o in enumerate(orig_ag, 1) if o in iface]
    repack(pose_ab, scorefxn, sub_ab)
    repack(pose_ag, scorefxn, sub_ag)
    e_ab, e_ag = scorefxn(pose_ab), scorefxn(pose_ag)

    jumps = vector1_int()
    setup_foldtree(pose_ia, spec, jumps)   # 重建 FoldTree，jump 1 分开两组链

    ia = InterfaceAnalyzerMover(1)
    ia.set_scorefunction(scorefxn)
    ia.set_pack_input(True)                # 两侧都 repack 才对称，缺一侧 dG 会偏移几十 REU
    ia.set_pack_separated(True)
    ia.set_calc_dSASA(True)
    ia.apply(pose_ia)
    t_ia = time.time() - t0

    return {
        'pdb_id': os.path.basename(path),
        'nres': pose_ia.total_residue(),
        'nres_repack': len(idx),
        'nres_antibody': pose_ab.total_residue(),
        'nres_antigen': pose_ag.total_residue(),
        'E_complex(REU)': round(e_complex, 2),
        'E_antibody(REU)': round(e_ab, 2),
        'E_antigen(REU)': round(e_ag, 2),
        'dG(REU)': round(e_complex - e_ab - e_ag, 2),      # 自己算，三者必然自洽
        'dG_IA(REU)': round(ia.get_interface_dG(), 2),     # InterfaceAnalyzer 的值，供对照
        'dSASA_int(A^2)': round(ia.get_interface_delta_sasa(), 1),
        'repack_time(s)': round(t_repack),
        'analyze_time(s)': round(t_ia),
    }


def done_set(out_csv):
    """已算过的文件名，用于断点续跑"""
    if not os.path.exists(out_csv):
        return set()
    with open(out_csv, newline='', encoding='utf-8') as f:
        return {row['pdb_id'] for row in csv.DictReader(f)}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--inputs', default=INPUTS)
    ap.add_argument('--out', default=OUTPUT)
    ap.add_argument('--interface', default=INTERFACE, help='链分组，如 HL_A')
    ap.add_argument('--radius', type=float, default=RADIUS)
    args = ap.parse_args()

    pyrosetta.init('-mute all')
    scorefxn = pyrosetta.get_fa_scorefxn()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    pdbs = sorted(glob.glob(os.path.join(args.inputs, '*.pdb')))
    skip = done_set(args.out)
    todo = [p for p in pdbs if os.path.basename(p) not in skip]

    print(f'界面 {args.interface}  |  共 {len(pdbs)} 个 pdb，'
          f'已完成 {len(skip)} 个，本次处理 {len(todo)} 个')

    new_file = not os.path.exists(args.out)
    with open(args.out, 'a', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        if new_file:
            w.writeheader()

        for n, path in enumerate(todo, 1):
            name = os.path.basename(path)
            try:
                row = evaluate(path, args.interface, args.radius, scorefxn)
                w.writerow(row)
                f.flush()        # 每条都落盘，中断不丢已完成的
                print(f"[{n}/{len(todo)}] {name}  dG={row['dG(REU)']}  "
                      f"({row['repack_time(s)'] + row['analyze_time(s)']} s)")
            except Exception as e:
                print(f'[{n}/{len(todo)}] {name}  失败: {type(e).__name__}: {e}')

    print('输出:', args.out)


if __name__ == '__main__':
    main()
