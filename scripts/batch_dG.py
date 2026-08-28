"""批量计算复合物的结合能 ΔG —— 不做任何 relax 与 repack。

直接对输入结构打分 -> 按链拆开 -> 直接打分 -> 相减
作为 batch_dG_repack_relax.py 的对照基线：数值不可靠，但快 50 倍
"""

# ==================== 配置（日常改这里） ====================

INPUTS    = '/data/lmk/rosetta_inputs'                          # 待评估的复合物 pdb 目录
OUTPUT    = '/data/lmk/rosetta_outputs/dG_results.csv'  # 指标输出
INTERFACE = 'HL_A'      # 链分组，要与 pdb 实际链号一致
                        # ⚠️ 左边必须是抗体、右边是抗原，否则 CSV 里两侧的列名会对调
# 01 不做任何 repack，这两套判据只影响报告列与 IA 的统计口径，详见 interface_lib.py
SITE_CRIT = 'heavy'     # 报告「真实界面有多大」—— 重原子接触
SITE_DIST = 4.0
PACK_DIST = 8.0         # 传给 InterfaceAnalyzer，与 02/03 的 repack 范围口径一致

# ===========================================================

import argparse
import csv
import glob
import os
import time

import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pyrosetta
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.protocols.docking import setup_foldtree
from pyrosetta.rosetta.utility import vector1_int

from interface_lib import group_pose, interface_both

FIELDS = ['pdb_id', 'residues_num_total', 'residues_num_interface',
          'residues_num_antibody', 'residues_num_antigen',
          'E_complex(REU)', 'E_antibody(REU)', 'E_antigen(REU)',
          'dG(REU)', 'dG_IA(REU)', 'dSASA_int(A^2)', 'total_time(s)']


def evaluate(path, spec, cfg, scorefxn):
    """算一个结构，返回一行指标。结构原样使用，不做任何优化"""
    t0 = time.time()
    pose = pyrosetta.pose_from_pdb(path)

    site = interface_both(pose, spec, cfg['site_dist'], cfg['site_crit'])   # 仅用于报告
    g1, g2 = spec.split('_')
    pose_ab, _ = group_pose(pose, g1)          # 直接拆开，不做任何优化
    pose_ag, _ = group_pose(pose, g2)

    e_complex = scorefxn(pose)
    e_ab, e_ag = scorefxn(pose_ab), scorefxn(pose_ag)

    jumps = vector1_int()
    setup_foldtree(pose, spec, jumps)    # 重建 FoldTree，jump 1 分开两组链

    ia = InterfaceAnalyzerMover(1)       # 独立走一遍，给出 dG_IA 与 dSASA
    ia.set_scorefunction(scorefxn)
    ia.set_pack_input(False)             # 复合物不 repack
    ia.set_pack_separated(False)         # 分开后也不 repack
    ia.set_calc_dSASA(True)

    ia.apply(pose)

    return {
        'pdb_id': os.path.basename(path),
        'residues_num_total': pose.total_residue(),
        'residues_num_interface': len(site),     # 重原子判据，真实界面大小
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
    ap.add_argument('--pack-dist', type=float, default=PACK_DIST)
    args = ap.parse_args()

    cfg = {'site_crit': args.site_criterion, 'site_dist': args.site_dist}
    # 与 02 / 03 的 repack 判据口径一致，dSASA 才可比
    pyrosetta.init(f'-mute all -pose_metrics:interface_cutoff {args.pack_dist}')
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
