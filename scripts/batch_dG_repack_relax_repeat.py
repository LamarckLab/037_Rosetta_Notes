"""生产脚本：批量算结合能 ΔG —— 重复采样 + 多进程并行。

在 03 的基础上做两件事：
  1. 每个结构重复 NSTRUCT 遍，给出 ΔG 的均值与标准差，而不是单次的一个数
  2. 按结构切分到 NPROC 个进程上并行

ΔG 取 InterfaceAnalyzer 的值。03 之所以另外自算一份，是为了拿到抗体、抗原
各自的能量；这里是汇总表，不报那两列，所以只留官方实现这一条路。

为什么要重复采样：FastRelax 每次落点不同，单次的 ΔG 不足以给设计排序 ——
实测同一结构两次 relax 的 E_complex 能差 100 REU 以上。02 的 repack 收敛得
很紧（干净结构逐位相同），所以只包 03。
"""

# ==================== 配置（日常改这里） ====================

INPUTS    = '/data/lmk/rosetta_inputs'                              # 待评估的复合物 pdb 目录
OUTPUT    = '/data/lmk/rosetta_outputs/dG_repack_relax_repeat_results.csv'     # 指标输出
INTERFACE = 'HL_A'      # 链分组，要与 pdb 实际链号一致
                        # ⚠️ 左边必须是抗体、右边是抗原，否则 CSV 里两侧的列名会对调
# 两套判据，回答两个不同的问题，详见 interface_lib.py
SITE_CRIT = 'heavy'     # 报告「真实界面有多大」—— 接触问题
SITE_DIST = 4.0
PACK_CRIT = 'cb'        # 决定 relax 与 repack 范围 —— 堆积问题
PACK_DIST = 8.0
NSTRUCT   = 5           # 每个结构重复几遍
NPROC     = 32          # 最多同时跑几个进程；实际取 min(NPROC, 待算结构数)

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

from batch_dG_repack_relax import evaluate    # 采样一遍 = 完整跑一遍 03，流程必须一致

FIELDS = ['pdb_id', 'residues_num_total', 'residues_num_interface',
          'residues_num_repack', 'residues_num_antibody', 'residues_num_antigen',
          'nstruct',
          'dG_mean(REU)', 'dG_std(REU)', 'dG_min(REU)', 'dG_max(REU)',
          'dG_all(REU)', 'dSASA_int_mean(A^2)', 'total_time(s)']

_SCOREFXN = None    # 每个 worker 自己的 scorefxn，不能跨进程传


def sample(path, spec, cfg, scorefxn, nstruct):
    """把 03 的 evaluate 跑 nstruct 遍，汇总成一行"""
    t0 = time.time()
    runs = [evaluate(path, spec, cfg, scorefxn) for _ in range(nstruct)]
    dgs = [r['dG_IA(REU)'] for r in runs]    # 取 IA 的值，不用 03 自算的那条
    first = runs[0]

    return {
        'pdb_id': first['pdb_id'],
        'residues_num_total': first['residues_num_total'],
        'residues_num_interface': first['residues_num_interface'],
        'residues_num_repack': first['residues_num_repack'],
        'residues_num_antibody': first['residues_num_antibody'],
        'residues_num_antigen': first['residues_num_antigen'],
        'nstruct': nstruct,
        'dG_mean(REU)': round(statistics.fmean(dgs), 2),
        'dG_std(REU)': round(statistics.stdev(dgs), 2) if nstruct > 1 else 0.0,
        'dG_min(REU)': min(dgs),
        'dG_max(REU)': max(dgs),
        'dG_all(REU)': ';'.join(str(v) for v in dgs),    # 原始值，留着重新统计
        'dSASA_int_mean(A^2)': round(
            statistics.fmean(r['dSASA_int(A^2)'] for r in runs), 1),
        'total_time(s)': round(time.time() - t0, 1),
    }


def init_worker(pack_dist):
    """每个 worker 起来时跑一次：独立 init，建自己的 scorefxn"""
    global _SCOREFXN
    import pyrosetta
    # silent=True 压掉 PyRosetta 的启动 banner，否则每个 worker 都刷一遍
    pyrosetta.init(f'-mute all -pose_metrics:interface_cutoff {pack_dist}', silent=True)
    _SCOREFXN = pyrosetta.get_fa_scorefxn()
    seed = pyrosetta.rosetta.numeric.random.rg().get_seed()
    print(f'  worker pid={os.getpid()}  seed={seed}', flush=True)


def work(task):
    """算一个结构。异常不能往外抛，否则整个 pool 会挂；带上 traceback 才好排查"""
    path, spec, cfg, nstruct = task
    try:
        return sample(path, spec, cfg, _SCOREFXN, nstruct)
    except Exception:
        return {'pdb_id': os.path.basename(path), '_error': traceback.format_exc()}


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
    ap.add_argument('--nstruct', type=int, default=NSTRUCT, help='每个结构重复几遍')
    ap.add_argument('--nproc', type=int, default=NPROC, help='最多同时跑几个进程')
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    pdbs = sorted(glob.glob(os.path.join(args.inputs, '*.pdb')))
    skip = done_set(args.out)
    todo = [p for p in pdbs if os.path.basename(p) not in skip]
    if not todo:
        print('没有待算的结构')
        return

    nproc = max(1, min(args.nproc, len(todo)))
    print(f'界面 {args.interface}  |  每个结构跑 {args.nstruct} 遍  |  {nproc} 进程并行\n'
          f'共 {len(pdbs)} 个 pdb，已完成 {len(skip)} 个，本次处理 {len(todo)} 个')

    cfg = {'site_crit': args.site_criterion, 'site_dist': args.site_dist,
           'pack_crit': args.pack_criterion, 'pack_dist': args.pack_dist}
    tasks = [(p, args.interface, cfg, args.nstruct) for p in todo]
    t0 = time.time()
    # 判「存在」不够：0 字节的残留文件会让表头永远写不出去，之后续跑会把首行数据读成字段名
    new_file = not (os.path.exists(args.out) and os.path.getsize(args.out) > 0)

    # spawn 而不是 fork：worker 拿到干净的解释器，不继承父进程里任何 Rosetta 状态
    ctx = mp.get_context('spawn')
    with open(args.out, 'a', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        if new_file:
            w.writeheader()

        with ctx.Pool(nproc, initializer=init_worker,
                      initargs=(args.pack_dist,)) as pool:
            for n, row in enumerate(pool.imap_unordered(work, tasks), 1):
                if '_error' in row:
                    print(f"[{n}/{len(todo)}] {row['pdb_id']}  失败:\n{row['_error']}",
                          flush=True)
                    continue
                w.writerow(row)
                f.flush()        # 每条都落盘，中断不丢已完成的
                print(f"[{n}/{len(todo)}] {row['pdb_id']}  "
                      f"dG={row['dG_mean(REU)']}±{row['dG_std(REU)']}  "
                      f"({row['total_time(s)']} s)", flush=True)

    print(f'\n总墙钟 {time.time() - t0:.1f} s')
    print('输出:', args.out)


if __name__ == '__main__':
    main()
