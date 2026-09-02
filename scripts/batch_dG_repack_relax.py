"""批量计算复合物的结合能 ΔG —— 界面 relax + 分开后 repack。

选界面残基 -> 只 relax 这些侧链 -> 打分 -> 拆开两侧各自 repack -> 打分 -> 相减
"""

# ==================== 配置（日常改这里） ====================

INPUTS    = '/data/lmk/rosetta_inputs'                          # 待评估的复合物 pdb 目录
OUTPUT    = '/data/lmk/rosetta_outputs/dG_repack_relax_results.csv'    # 指标输出
INTERFACE = 'HL_A'      # 链分组，要与 pdb 实际链号一致
                        # ⚠️ 左边必须是抗体、右边是抗原，否则 CSV 里两侧的列名会对调
# 两套判据，回答两个不同的问题，详见 interface_lib.py
SITE_CRIT = 'heavy'     # 报告「真实界面有多大」—— 接触问题
SITE_DIST = 4.0
PACK_CRIT = 'cb'        # 决定 relax 与 repack 范围 —— 堆积问题
PACK_DIST = 8.0         # 同时会传给 InterfaceAnalyzer，两条路的 repack 范围才可比

# ===========================================================

import argparse
import csv
import glob
import os
import time

import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pyrosetta
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.protocols.docking import setup_foldtree
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.utility import vector1_int

from interface_lib import group_pose, interface_both, repack

FIELDS = ['pdb_id', 'residues_num_total', 'residues_num_interface',
          'residues_num_repack', 'residues_num_antibody', 'residues_num_antigen',
          'E_complex(REU)', 'E_antibody(REU)', 'E_antigen(REU)',
          'dG(REU)', 'dG_IA(REU)', 'dSASA_int(A^2)', 'total_time(s)']




def evaluate(path, spec, cfg, scorefxn):
    """算一个结构，返回一行指标。并行时每个进程各自 init 并传入自己的 scorefxn"""
    t0 = time.time()
    pose = pyrosetta.pose_from_pdb(path)

    idx = interface_both(pose, spec, cfg['pack_dist'], cfg['pack_crit'])   # relax+repack 范围
    site = interface_both(pose, spec, cfg['site_dist'], cfg['site_crit'])  # 仅用于报告
    if not idx:
        raise ValueError(f'界面为空，检查链分组 {spec} 是否与该文件匹配')

    mm = MoveMap()              # 主链固定，只放开界面残基的侧链
    mm.set_bb(False)
    mm.set_chi(False)
    for i in idx:
        mm.set_chi(i, True)

    fr = FastRelax()
    fr.set_scorefxn(scorefxn)
    fr.set_movemap(mm)

    fr.apply(pose)

    pose_ia = pose.clone()           # 留一份 relax 后的原样，供 IA 独立走一遍完整流程

    repack(pose, scorefxn, idx)      # 结合态也做一次同范围 repack，与分离态对等
    e_complex = scorefxn(pose)

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

    ia = InterfaceAnalyzerMover(1)         # 独立走一遍，给出 dG_IA 与 dSASA
    ia.set_scorefunction(scorefxn)
    ia.set_pack_input(True)                # 两侧都 repack 才对称，缺一侧 dG 会偏移几十 REU
    ia.set_pack_separated(True)
    ia.set_calc_dSASA(True)
    ia.apply(pose_ia)

    return {
        'pdb_id': os.path.basename(path),
        'residues_num_total': pose_ia.total_residue(),
        'residues_num_interface': len(site),     # 重原子判据，真实界面大小
        'residues_num_repack': len(idx),         # CB 判据，实际放开的侧链
        'residues_num_antibody': pose_ab.total_residue(),
        'residues_num_antigen': pose_ag.total_residue(),
        'E_complex(REU)': round(e_complex, 2),
        'E_antibody(REU)': round(e_ab, 2),
        'E_antigen(REU)': round(e_ag, 2),
        'dG(REU)': round(e_complex - e_ab - e_ag, 2),      # 与前三列必然自洽，可直接验算
        'dG_IA(REU)': round(ia.get_interface_dG(), 2),     # InterfaceAnalyzer 的值，供对照
        'dSASA_int(A^2)': round(ia.get_interface_delta_sasa(), 1),
        'total_time(s)': round(time.time() - t0, 1),
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
    ap.add_argument('--site-criterion', default=SITE_CRIT, choices=['cb', 'heavy'])
    ap.add_argument('--site-dist', type=float, default=SITE_DIST)
    ap.add_argument('--pack-criterion', default=PACK_CRIT, choices=['cb', 'heavy'])
    ap.add_argument('--pack-dist', type=float, default=PACK_DIST)
    args = ap.parse_args()

    cfg = {'site_crit': args.site_criterion, 'site_dist': args.site_dist,
           'pack_crit': args.pack_criterion, 'pack_dist': args.pack_dist}
    # IA 的 repack 范围也设成同一个值，否则 dG 与 dG_IA 又会因判据不同而分歧
    pyrosetta.init(f'-mute all -pose_metrics:interface_cutoff {args.pack_dist}')
    scorefxn = pyrosetta.get_fa_scorefxn()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    pdbs = sorted(glob.glob(os.path.join(args.inputs, '*.pdb')))
    skip = done_set(args.out)
    todo = [p for p in pdbs if os.path.basename(p) not in skip]

    print(f'界面 {args.interface}  |  共 {len(pdbs)} 个 pdb，'
          f'已完成 {len(skip)} 个，本次处理 {len(todo)} 个')

    # 判「存在」不够：0 字节的残留文件会让表头永远写不出去，之后续跑会把首行数据读成字段名
    new_file = not (os.path.exists(args.out) and os.path.getsize(args.out) > 0)
    with open(args.out, 'a', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        if new_file:
            w.writeheader()

        for n, path in enumerate(todo, 1):
            name = os.path.basename(path)
            try:
                row = evaluate(path, args.interface, cfg, scorefxn)
                w.writerow(row)
                f.flush()        # 每条都落盘，中断不丢已完成的
                print(f"[{n}/{len(todo)}] {name}  dG={row['dG(REU)']}  "
                      f"({row['total_time(s)']} s)")
            except Exception as e:
                print(f'[{n}/{len(todo)}] {name}  失败: {type(e).__name__}: {e}')

    print('输出:', args.out)


if __name__ == '__main__':
    main()
