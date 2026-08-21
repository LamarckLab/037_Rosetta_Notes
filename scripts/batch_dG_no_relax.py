"""批量计算复合物的结合能 ΔG —— 不做任何 relax 与 repack。

直接对输入结构打分 -> 沿 jump 拆开 -> 直接打分 -> 相减
作为 batch_dG_relax.py 的对照基线：数值不可靠，但快 50 倍
"""

# ==================== 配置（日常改这里） ====================

INPUTS    = '/data/lmk/rosetta_inputs'                          # 待评估的复合物 pdb 目录
OUTPUT    = '/data/lmk/rosetta_outputs/dG_no_relax_results.csv'  # 指标输出
INTERFACE = 'HL_A'      # 链分组，要与 pdb 实际链号一致
                        # ⚠️ 左边必须是抗体、右边是抗原，否则 CSV 里两侧的列名会对调

# ===========================================================

import argparse
import csv
import glob
import os
import time

import pyrosetta
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.protocols.docking import setup_foldtree
from pyrosetta.rosetta.utility import vector1_int

FIELDS = ['pdb_id', 'nres', 'nres_antibody', 'nres_antigen',
          'E_complex(REU)', 'E_antibody(REU)', 'E_antigen(REU)', 'dG(REU)',
          'dSASA_int(A^2)', 'analyze_time(s)']


def evaluate(path, spec, scorefxn):
    """算一个结构，返回一行指标。结构原样使用，不做任何优化"""
    pose = pyrosetta.pose_from_pdb(path)

    jumps = vector1_int()
    setup_foldtree(pose, spec, jumps)    # 重建 FoldTree，jump 1 分开两组链

    ia = InterfaceAnalyzerMover(1)
    ia.set_scorefunction(scorefxn)
    ia.set_pack_input(False)             # 复合物不 repack
    ia.set_pack_separated(False)         # 分开后也不 repack
    ia.set_calc_dSASA(True)

    t0 = time.time()
    ia.apply(pose)
    t_ia = time.time() - t0

    return {
        'pdb_id': os.path.basename(path),
        'nres': pose.total_residue(),
        'nres_antibody': ia.get_side1_nres(),       # side1 = INTERFACE 下划线左边那组
        'nres_antigen': ia.get_side2_nres(),        # side2 = 右边那组
        'E_complex(REU)': round(ia.get_complex_energy(), 2),
        'E_antibody(REU)': round(ia.get_side1_score(), 2),
        'E_antigen(REU)': round(ia.get_side2_score(), 2),
        'dG(REU)': round(ia.get_interface_dG(), 2),
        'dSASA_int(A^2)': round(ia.get_interface_delta_sasa(), 1),
        'analyze_time(s)': round(t_ia, 1),
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
                row = evaluate(path, args.interface, scorefxn)
                w.writerow(row)
                f.flush()        # 每条都落盘，中断不丢已完成的
                print(f"[{n}/{len(todo)}] {name}  dG={row['dG(REU)']}  "
                      f"({row['analyze_time(s)']} s)")
            except Exception as e:
                print(f'[{n}/{len(todo)}] {name}  失败: {type(e).__name__}: {e}')

    print('输出:', args.out)


if __name__ == '__main__':
    main()
