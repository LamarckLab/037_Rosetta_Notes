"""第一级筛选：对界面抗体侧做饱和突变，快速排出 ΔΔG 的高低顺序。

对每个界面位置枚举其余 19 种氨基酸，各按 02 的流程（REF2015，只 repack）算一次
ΔG，与野生型相减。

⚠️ 这一列叫 ddG_screen 而不是 ddG，是有意的：REF2015 + 只 repack 不是界面 ΔΔG
的标准算法，突变成更大的残基时无法靠主链微调化解冲突，会被系统性高估。**它只能
用来排序，不能作为 ΔΔG 报告。** 报告值来自第二级的 flex ddG（talaris2014）。
两级的数字隔着一整个能量函数，任何情况下都不要放在一起比。

符号：ddG_screen < 0 表示突变体结合更强。
"""

# ==================== 配置（日常改这里） ====================

INPUTS    = '/data/lmk/rosetta_inputs'                          # 待扫描的复合物 pdb 目录
OUTPUT    = '/data/lmk/rosetta_outputs/ddG_screen_results.csv'  # 指标输出
INTERFACE = 'HL_A'      # 链分组，要与 pdb 实际链号一致
                        # ⚠️ 左边必须是抗体、右边是抗原
RADIUS    = 10.0        # 界面判定半径 A（CB 之间的距离）
MUT_SIDE  = 'antibody'  # 突变哪一侧：antibody（下划线左边）或 antigen
WT_NREPEAT= 20          # 野生型参考值重复几次
WT_TRIM   = 5           # 排序后掐掉最负的 N 个和最正的 N 个，用中间的求均值
                        # packer 约 20% 概率掉进一个差 3.5 REU 的次优解，而
                        # dG_wt 锚定全部突变，单次抽样会把整张表平移
NPROC     = 32          # 最多同时跑几个进程

# ===========================================================

import argparse
import csv
import glob
import multiprocessing as mp
import os
import statistics
import sys
import time
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from batch_dG_repack import (chains_selector, group_pose, interface_residues, repack)

AA3 = {'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
       'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
       'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
       'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'}

FIELDS = ['pdb_id', 'chain', 'pdb_position', 'icode', 'wt_aa', 'mut_aa',
          'dG_wt(REU)', 'dG_mut(REU)', 'ddG_screen(REU)', 'rank']

_SCOREFXN = None    # 每个 worker 自己的 scorefxn，不能跨进程传
_POSE = {}          # 每个 worker 缓存已读过的 pdb，避免每个突变都重读


def mutate_side(pose, spec, radius, side):
    """要突变的那一侧的界面残基编号（Rosetta 编号）"""
    from pyrosetta.rosetta.core.select import get_residues_from_subset
    from pyrosetta.rosetta.core.select.residue_selector import (
        AndResidueSelector, NeighborhoodResidueSelector)

    g1, g2 = spec.split('_')
    mine, other = (g1, g2) if side == 'antibody' else (g2, g1)
    sel = AndResidueSelector(
        NeighborhoodResidueSelector(chains_selector(other), radius, False),
        chains_selector(mine))
    return sorted(get_residues_from_subset(sel.apply(pose)))


def dG_from_pose(pose, spec, idx, scorefxn):
    """02 的流程，但入参是内存里的 Pose 而不是文件路径"""
    repack(pose, scorefxn, idx)                  # 结合态 repack
    e_complex = scorefxn(pose)

    g1, g2 = spec.split('_')
    iface = set(idx)
    pose_ab, orig_ab = group_pose(pose, g1)      # 拆开
    pose_ag, orig_ag = group_pose(pose, g2)
    repack(pose_ab, scorefxn,
           [j for j, o in enumerate(orig_ab, 1) if o in iface])
    repack(pose_ag, scorefxn,
           [j for j, o in enumerate(orig_ag, 1) if o in iface])
    return e_complex - scorefxn(pose_ab) - scorefxn(pose_ag)


def load(path):
    """读 pdb 并缓存；同一个 worker 处理同一结构的上百个突变，只读一次"""
    import pyrosetta
    if path not in _POSE:
        _POSE[path] = pyrosetta.pose_from_pdb(path)
    return _POSE[path]


def wt_setup(path, spec, radius, side):
    """界面残基与待突变位点清单，不算能量"""
    pose = load(path)
    idx = interface_residues(pose, spec, radius)         # repack 范围：界面两侧
    sites = mutate_side(pose, spec, radius, side)        # 突变范围：只一侧
    if not sites:
        raise ValueError(f'界面为空，检查链分组 {spec} 与 MUT_SIDE={side}')

    info = pose.pdb_info()
    return {
        'repack_idx': idx,
        # icode 必须带上：H100 与 H100A 是两个不同的残基，只记 number 会重名
        'sites': [(p, info.chain(p), info.number(p), info.icode(p).strip(),
                   pose.residue(p).name1()) for p in sites],
    }


def wt_energy(task):
    """野生型 ΔG 的一次抽样，重复多次由主进程汇总"""
    path, spec, idx = task
    return dG_from_pose(load(path).clone(), spec, idx, _SCOREFXN)


def trimmed(vals, trim):
    """排序后掐掉两端各 trim 个，返回中间部分的均值"""
    v = sorted(vals)
    kept = v[trim:len(v) - trim] or v
    return round(statistics.fmean(kept), 2), kept


def init_worker(radius):
    """每个 worker 起来时跑一次：独立 init，建自己的 scorefxn"""
    global _SCOREFXN
    import pyrosetta
    pyrosetta.init(f'-mute all -pose_metrics:interface_cutoff {radius}', silent=True)
    _SCOREFXN = pyrosetta.get_fa_scorefxn()


def work(task):
    """算一个突变。异常不能往外抛，否则整个 pool 会挂"""
    path, spec, pos, chain, num, icode, wt_aa, mut_aa, idx, dG_wt = task
    try:
        from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
        pose = load(path).clone()
        MutateResidue(pos, AA3[mut_aa]).apply(pose)
        dG_mut = dG_from_pose(pose, spec, idx, _SCOREFXN)
        return {
            'pdb_id': os.path.basename(path),
            'chain': chain, 'pdb_position': num, 'icode': icode,
            'wt_aa': wt_aa, 'mut_aa': mut_aa,
            'dG_wt(REU)': dG_wt,
            'dG_mut(REU)': round(dG_mut, 2),
            'ddG_screen(REU)': round(dG_mut - dG_wt, 2),
            'rank': '',                      # 阶段三统一填
        }
    except Exception:
        return {'pdb_id': os.path.basename(path),
                '_error': f'{chain}{num}{icode}{wt_aa}->{mut_aa}\n{traceback.format_exc()}'}


def key(row):
    """断点续跑用的唯一键"""
    return (f"{row['pdb_id']}:{row['chain']}{row['pdb_position']}"
            f"{row.get('icode', '')}{row['mut_aa']}")


def done_set(out_csv):
    if not os.path.exists(out_csv):
        return set()
    with open(out_csv, newline='', encoding='utf-8') as f:
        return {key(r) for r in csv.DictReader(f)}


def rerank(out_csv):
    """全部算完后按 ddG_screen 升序重排，负得越多排越前"""
    with open(out_csv, newline='', encoding='utf-8') as f:
        rows = list(csv.DictReader(f))
    rows.sort(key=lambda r: float(r['ddG_screen(REU)']))
    for n, r in enumerate(rows, 1):
        r['rank'] = n
    with open(out_csv, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows({k: r[k] for k in FIELDS} for r in rows)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--inputs', default=INPUTS)
    ap.add_argument('--out', default=OUTPUT)
    ap.add_argument('--interface', default=INTERFACE, help='链分组，如 HL_A')
    ap.add_argument('--radius', type=float, default=RADIUS)
    ap.add_argument('--side', default=MUT_SIDE, choices=['antibody', 'antigen'])
    ap.add_argument('--wt-nrepeat', type=int, default=WT_NREPEAT)
    ap.add_argument('--wt-trim', type=int, default=WT_TRIM)
    ap.add_argument('--nproc', type=int, default=NPROC)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    pdbs = sorted(glob.glob(os.path.join(args.inputs, '*.pdb')))
    skip = done_set(args.out)
    print(f'界面 {args.interface}  |  突变 {args.side} 侧  |  '
          f'{args.nproc} 进程  |  {len(pdbs)} 个 pdb，已完成 {len(skip)} 个突变')

    t0 = time.time()
    new_file = not os.path.exists(args.out)
    ctx = mp.get_context('spawn')

    with open(args.out, 'a', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        if new_file:
            w.writeheader()

        with ctx.Pool(args.nproc, initializer=init_worker,
                      initargs=(args.radius,)) as pool:
            for path in pdbs:
                name = os.path.basename(path)
                ref = pool.apply(wt_setup,
                                 (path, args.interface, args.radius, args.side))

                # 野生型锚定全部突变，多抽几次截尾平均，免得被一次坏解整体平移
                draws = pool.map(wt_energy,
                                 [(path, args.interface, ref['repack_idx'])]
                                 * args.wt_nrepeat)
                dG_wt, kept = trimmed(draws, args.wt_trim)
                print('')
                print(f'{name}  dG_wt={dG_wt}   {args.wt_nrepeat} 次抽样 '
                      f'[{min(draws):.2f}, {max(draws):.2f}] -> 保留 {len(kept)} 个 '
                      f'[{min(kept):.2f}, {max(kept):.2f}]')

                tasks = [(path, args.interface, pos, ch, num, ic, wt, mut,
                          ref['repack_idx'], dG_wt)
                         for pos, ch, num, ic, wt in ref['sites']
                         for mut in AA3 if mut != wt]
                tasks = [t for t in tasks
                         if f'{name}:{t[3]}{t[4]}{t[5]}{t[7]}' not in skip]
                print(f'{len(ref["sites"])} 个位点，{len(tasks)} 个待算突变')
                if not tasks:
                    continue

                for n, row in enumerate(pool.imap_unordered(work, tasks), 1):
                    if '_error' in row:
                        print(f"  失败 {row['_error']}", flush=True)
                        continue
                    w.writerow(row)
                    f.flush()        # 每条都落盘，中断不丢已完成的
                    if n % 50 == 0 or n == len(tasks):
                        print(f'  [{n}/{len(tasks)}]', flush=True)

    rows = rerank(args.out)
    print(f'\n总墙钟 {time.time() - t0:.1f} s，共 {len(rows)} 个突变')
    print('输出:', args.out)
    print('\n--- ddG_screen 最负的 10 个（仅供排序，不是 ΔΔG 报告值）---')
    for r in rows[:10]:
        print(f"  {r['pdb_id']:<16} "
              f"{r['chain']}{r['pdb_position']}{r['icode']:<5} "
              f"{r['wt_aa']} -> {r['mut_aa']}   {r['ddG_screen(REU)']:>8}")


if __name__ == '__main__':
    main()
