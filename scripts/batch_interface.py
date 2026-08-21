"""批量计算复合物的界面结合能 dG_separated。

读 PDB -> 只 relax 界面残基的侧链 -> InterfaceAnalyzer -> 追加一行到 CSV
"""

# ==================== 配置（日常改这里） ====================

INPUTS    = '/data/lmk/rosetta_inputs'                          # 待评估的复合物 pdb 目录
OUTPUT    = '/data/lmk/rosetta_outputs/interface_metrics.csv'   # 指标输出
INTERFACE = 'HL_A'      # 链分组，下划线左右各一组；要与 pdb 里的实际链号一致
RADIUS    = 8.0         # 界面判定半径 A

# ===========================================================

import argparse
import csv
import glob
import os
import time

import pyrosetta
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.select import get_residues_from_subset
from pyrosetta.rosetta.core.select.residue_selector import (
    AndResidueSelector, ChainSelector, NeighborhoodResidueSelector, OrResidueSelector)
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.protocols.docking import setup_foldtree
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.utility import vector1_int

FIELDS = ['pdb', 'nres', 'nres_iface', 'dG_separated', 'dSASA',
          'total_score', 'relax_s', 'ia_s']


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


def evaluate(path, spec, radius, scorefxn):
    """算一个结构，返回一行指标。并行时每个进程各自 init 并传入自己的 scorefxn"""
    pose = pyrosetta.pose_from_pdb(path)

    idx = interface_residues(pose, spec, radius)
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

    t0 = time.time()
    fr.apply(pose)
    t_relax = time.time() - t0

    jumps = vector1_int()
    setup_foldtree(pose, spec, jumps)    # 重建 FoldTree，jump 1 分开两组链

    ia = InterfaceAnalyzerMover(1)
    ia.set_scorefunction(scorefxn)
    ia.set_pack_separated(True)          # 分开后重新 repack，模拟解离
    ia.set_calc_dSASA(True)              # dSASA 默认不算，要显式打开

    t0 = time.time()
    ia.apply(pose)
    t_ia = time.time() - t0

    return {
        'pdb': os.path.basename(path),
        'nres': pose.total_residue(),
        'nres_iface': len(idx),
        'dG_separated': round(ia.get_interface_dG(), 2),
        'dSASA': round(ia.get_interface_delta_sasa(), 1),
        'total_score': round(scorefxn(pose), 2),
        'relax_s': round(t_relax),
        'ia_s': round(t_ia),
    }


def done_set(out_csv):
    """已算过的文件名，用于断点续跑"""
    if not os.path.exists(out_csv):
        return set()
    with open(out_csv, newline='', encoding='utf-8') as f:
        return {row['pdb'] for row in csv.DictReader(f)}


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
                print(f"[{n}/{len(todo)}] {name}  dG={row['dG_separated']}  "
                      f"({row['relax_s'] + row['ia_s']} s)")
            except Exception as e:
                print(f'[{n}/{len(todo)}] {name}  失败: {type(e).__name__}: {e}')

    print('输出:', args.out)


if __name__ == '__main__':
    main()
