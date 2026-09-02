"""第二级的实现：对筛选后的突变跑官方 flex ddG，产出可报告的 ΔΔG。

不直接跑，由 batch_ddG_flex_antibody.py / _antigen.py 两个入口调用，它们只覆盖
默认的输入输出路径。

读第一级输出（人工筛过的）csv，逐个突变按 Kortemme 实验室的 ddG-backrub.xml
跑 NSTRUCT 条 backrub 轨迹，取数与重加权照官方 analyze_flex_ddG.py。

⚠️ 本档用 **talaris2014**，不是前几档的 REF2015 —— flex ddG 整套是在 talaris2014
上标定的，实测换 ref2015 相关性会掉（r 0.79 -> 0.57~0.68，Sora 等 2023）。这需要
整个进程开 -corrections::restore_talaris_behavior，所以只能独立进程跑。
**本档的数字与 01-04 的 REU 不可比。**

符号：ddG < 0 表示突变体结合更强。

主指标 `ddG_mean(REU)` 是七个分子间能量项双差之后的直接加和，未经任何标定。
`ddG_gam_*` 是同一批数据过 ZEMu GAM 重加权的结果，单位 kcal/mol，供对照 ——
论文报的 r = 0.79 / MAE ≈ 1 kcal/mol 是 GAM 那一路的成绩，且 GAM 是逐项过
sigmoid 再求和，**两种口径的排序可能不同**，不是等比缩放。

参考 https://github.com/Kortemm  e-Lab/flex_ddG_tutorial
"""

# ==================== 配置（日常改这里） ====================

SELECTED  = '/data/lmk/rosetta_outputs/ddG_selected.csv'        # 人工筛过的突变清单
PDB_DIR   = '/data/lmk/rosetta_inputs'                          # 按 pdb_id 到这里找结构
OUTPUT    = '/data/lmk/rosetta_outputs/ddG_flex_results.csv'    # 指标输出
XML       = '/data/lmk/rosetta_scripts/ddG-backrub.xml'         # 官方协议，原样保存
WORK      = '/data/lmk/rosetta_work/flexddg'                    # 每条轨迹的临时 db3

MOVE_CHAIN     = 'A'      # 算解离态时把哪条（组）链移开，一般就是抗原
BUBBLE         = 'cb8'    # backrub 支点与 repack 的邻域判据，围绕突变位点画：
                          #   cb8    官方原版，邻居原子(CB) 8 Å —— r=0.79 就是它标定的
                          #   atom4  任意重原子 4 Å；实测两者选出的残基数相当
                          #          （12~13 vs 13~15），但构成不同
                          # ⚠️ 改成 atom4 就不再是官方协议，论文的精度数据不适用
NSTRUCT        = 35       # 每个突变跑几条轨迹；官方推荐 35
BACKRUB_TRIALS = 35000    # 每条轨迹的 backrub 步数；官方推荐 35000
MAX_MIN_ITER   = 5000     # 官方 benchmark 值
CONV_THRESH    = 1.0      # abs_score_convergence_thresh，官方值
NPROC          = 32       # 最多同时跑几个进程

# ===========================================================

import argparse
import csv
import os
import statistics
import sys
import time
import traceback
import multiprocessing as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from flexddg_lib import ddG_from_db3

KEYS = ['pdb_id', 'chain', 'pdb_position', 'icode', 'wt_aa', 'mut_aa']
FIELDS = KEYS + ['ddG_mean(REU)', 'ddG_std(REU)', 'ddG_all(REU)',
                 'ddG_gam_mean(kcal/mol)', 'ddG_gam_std(kcal/mol)',
                 'nstruct_done', 'backrub_trials', 'cpu_time(s)']

# 官方 XML 里定义 bubble 的那一行；变体只替换它，磁盘上的文件始终与上游逐字节一致
BUBBLE_CB8 = '<Neighborhood name="bubble" selector="resselector" distance="8.0"/>'
BUBBLE_ATOM4 = ('<CloseContact name="bubble" residue_selector="resselector" '
                'contact_threshold="4.0"/>')

_XML = None      # 每个 worker 读一次模板
_POSE = {}       # 每个 worker 缓存已读过的 pdb


def tag_of(row):
    return (f"{row['pdb_id']}:{row['chain']}{row['pdb_position']}"
            f"{row.get('icode', '')}{row['mut_aa']}")


def locate(pose, chain, num, icode):
    """按 链+编号+插入编号 找 Rosetta 内部编号

    抗体 Kabat 编号里 H100 与 H100A 是两个不同的残基，只按 number 定位会张冠李戴。
    """
    info = pose.pdb_info()
    for i in range(1, pose.total_residue() + 1):
        if (info.chain(i) == chain and info.number(i) == int(num)
                and info.icode(i).strip() == (icode or '').strip()):
            return i
    return None


def init_worker(xml_path, work, opts, bubble):
    """每个 worker 起来时跑一次。注意这里进的是 talaris 模式，全局生效"""
    global _XML
    import pyrosetta
    # -inout:dbms:database_name 必须给个值（否则 ReportToDB 构建时报 inactive option），
    # 但真正写哪个库由 XML 里的 database_name 属性决定，逐条轨迹各写各的
    pyrosetta.init(
        '-mute all -corrections::restore_talaris_behavior true '
        '-inout:dbms:mode sqlite3 '
        f'-inout:dbms:database_name {work}/unused_{os.getpid()}.db3 '
        '-in:file:fullatom -ignore_unrecognized_res '
        '-ignore_zero_occupancy false -ex1 -ex2', silent=True)
    _XML = open(xml_path).read()
    if bubble == 'atom4':
        assert BUBBLE_CB8 in _XML, 'XML 里找不到官方那行 Neighborhood，无法替换'
        _XML = _XML.replace(BUBBLE_CB8, BUBBLE_ATOM4)


def load(path):
    import pyrosetta
    if path not in _POSE:
        _POSE[path] = pyrosetta.pose_from_pdb(path)
    return _POSE[path]


def run_trajectory(task):
    """跑一条 backrub 轨迹，返回这一条的 ΔΔG"""
    from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects

    row, k, pdb_path, work, opts = task
    t0 = time.time()
    tag = f"{row['pdb_id'][:-4]}_{row['chain']}{row['pdb_position']}{row['mut_aa']}_{k}"
    d = os.path.join(work, tag)
    try:
        os.makedirs(d, exist_ok=True)
        resfile = os.path.join(d, 'mutate.resfile')
        # 动手前先核对：csv 说这里是 wt_aa，结构里也必须真的是它。编号体系对不上时
        # （插入编号丢失、换了 pdb 版本、手工改错 csv），这里当场报错，而不是安静地
        # 算出一个张冠李戴的 ΔΔG
        pose = load(pdb_path).clone()
        icode = (row.get('icode') or '').strip()
        target = locate(pose, row['chain'], row['pdb_position'], icode)
        if target is None:
            raise ValueError(
                f"结构里找不到 {row['chain']}{row['pdb_position']}{icode}")
        actual = pose.residue(target).name1()
        if actual != row['wt_aa']:
            raise ValueError(
                f"{row['chain']}{row['pdb_position']}{icode} 实际是 {actual}，"
                f"csv 写的是 {row['wt_aa']} —— 编号对不上，拒绝计算")

        # 官方 resfile 格式；插入编号直接跟在数字后面无空格（已实测 100A 能精确命中）
        mut_line = (f"{row['pdb_position']}{icode} {row['chain']} "
                    f"PIKAA {row['mut_aa']}")
        with open(resfile, 'w') as f:
            f.write(chr(10).join(['NATAA', 'start', mut_line, '']))

        db3 = os.path.join(d, 'ddG.db3')
        xml = (_XML.replace('database_name="ddG.db3"', f'database_name="{db3}"')
                   .replace('database_name="struct.db3"',
                            f'database_name="{os.path.join(d, "struct.db3")}"'))
        for key, val in {'mutate_resfile_relpath': resfile,
                         'chainstomove': opts['move_chain'],
                         'number_backrub_trials': str(opts['trials']),
                         'backrub_trajectory_stride': str(opts['trials']),
                         'max_minimization_iter': str(opts['min_iter']),
                         'abs_score_convergence_thresh': str(opts['conv'])}.items():
            xml = xml.replace(f'%%{key}%%', val)

        XmlObjects.create_from_string(xml).get_mover('ParsedProtocol').apply(pose)
        gam, raw, _ = ddG_from_db3(db3)
        res = {'key': tag_of(row), 'gam': gam, 'raw': raw,
               'sec': time.time() - t0}
    except Exception:
        return {'key': tag_of(row), '_error': traceback.format_exc(), 'dir': d}

    # 轨迹目录一律保留，成功失败都不删。单条约 1.6 MB，700 条也才 1 GB 出头，
    # 而每条是 31 分钟 CPU 换来的 —— 删掉就只剩 csv 里的聚合值，想换统计口径
    # 或事后排查都得重跑。失败时它更是唯一的现场。
    return res


def done_set(out_csv):
    if not os.path.exists(out_csv):
        return set()
    with open(out_csv, newline='', encoding='utf-8') as f:
        return {tag_of(r) for r in csv.DictReader(f)}


def summarize(row, vals, trials):
    """把一个突变的若干条轨迹汇总成一行

    vals 为空也要出一行 —— 一个突变一条轨迹都没取到数时（多半是 selected.csv
    里的位置在 pdb 里不存在），指标列留空、nstruct_done=0，这样输出行数与输入
    永远对得上，缺哪个一眼可见，不用回头翻日志。
    """
    base = {**{k: row[k] for k in KEYS},
            'nstruct_done': len(vals),
            'backrub_trials': trials,
            # 各条轨迹耗时之和，即该突变真实占用的 CPU 时间；并行下墙钟没有可比性
            'cpu_time(s)': round(sum(v['sec'] for v in vals), 1)}
    if not vals:
        return {**base, **{k: '' for k in FIELDS if k.startswith('ddG')}}

    gams = [v['gam'] for v in vals]
    raws = [v['raw'] for v in vals]
    two = len(vals) > 1
    return {
        **base,
        'ddG_mean(REU)': round(statistics.fmean(raws), 3),
        'ddG_std(REU)': round(statistics.stdev(raws), 3) if two else 0.0,
        # 每条轨迹的原始值，留着换统计口径（中位数、截尾均值、看分布）时不用重跑。
        # 一条轨迹是 31 分钟 CPU，重算 35 条要 18 小时，这一列必须留
        'ddG_all(REU)': ';'.join(f'{r:.3f}' for r in raws),
        'ddG_gam_mean(kcal/mol)': round(statistics.fmean(gams), 3),
        'ddG_gam_std(kcal/mol)': round(statistics.stdev(gams), 3) if two else 0.0,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--selected', default=SELECTED, help='人工筛过的突变 csv')
    ap.add_argument('--pdb-dir', default=PDB_DIR)
    ap.add_argument('--out', default=OUTPUT)
    ap.add_argument('--xml', default=XML)
    ap.add_argument('--work', default=WORK)
    ap.add_argument('--move-chain', default=MOVE_CHAIN, help='解离时移开哪条链')
    ap.add_argument('--bubble', default=BUBBLE, choices=['cb8', 'atom4'],
                    help='突变位点邻域判据；cb8 为官方原版')
    ap.add_argument('--nstruct', type=int, default=NSTRUCT)
    ap.add_argument('--trials', type=int, default=BACKRUB_TRIALS)
    ap.add_argument('--nproc', type=int, default=NPROC)
    args = ap.parse_args()

    os.makedirs(args.work, exist_ok=True)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    with open(args.selected, newline='', encoding='utf-8') as f:
        rows = [r for r in csv.DictReader(f)]
    skip = done_set(args.out)
    rows = [r for r in rows if tag_of(r) not in skip]
    if not rows:
        print('没有待算的突变')
        return

    opts = {'move_chain': args.move_chain, 'trials': args.trials,
            'min_iter': MAX_MIN_ITER, 'conv': CONV_THRESH}
    # 任务粒度 = 一条轨迹。若按突变切，每个任务几十小时，进程之间没法均衡
    tasks = [(r, k, os.path.join(args.pdb_dir, r['pdb_id']), args.work, opts)
             for r in rows for k in range(args.nstruct)]
    by_key = {tag_of(r): r for r in rows}

    print(f'{len(rows)} 个突变 × {args.nstruct} 条轨迹 = {len(tasks)} 个任务  |  '
          f'{args.trials} backrub steps  |  bubble {args.bubble}  |  {args.nproc} 进程\n'
          f'⚠️ talaris2014 模式，本档结果与 01-04 的 REU 不可比')

    t0 = time.time()
    got = {}                 # key -> 成功的轨迹结果
    fails = {}               # key -> [失败次数, 首条 traceback, 现场目录]
    written = set()          # 已经落盘的突变，收尾时据此补齐剩下的
    # 判「存在」不够：0 字节的残留文件会让表头永远写不出去，之后续跑会把首行数据读成字段名
    new_file = not (os.path.exists(args.out) and os.path.getsize(args.out) > 0)
    ctx = mp.get_context('spawn')

    with open(args.out, 'a', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        if new_file:
            w.writeheader()

        with ctx.Pool(args.nproc, initializer=init_worker,
                      initargs=(args.xml, args.work, opts, args.bubble)) as pool:
            for n, res in enumerate(pool.imap_unordered(run_trajectory, tasks), 1):
                key = res['key']
                if '_error' in res:
                    # 首条打完整 traceback，其余只累计 —— 但计数必须留着，
                    # 否则「一个突变全灭」在日志上会显示成「偶发一例」
                    if key not in fails:
                        fails[key] = [0, res['_error'], res.get('dir', '?')]
                        print(f"[{n}/{len(tasks)}] {key} 首次失败，"
                              f"现场保留在 {res.get('dir', '?')}", flush=True)
                        print(res['_error'], flush=True)
                    fails[key][0] += 1
                    continue

                got.setdefault(key, []).append(res)
                if len(got[key]) < args.nstruct:
                    continue

                row = summarize(by_key[key], got.pop(key), args.trials)
                written.add(key)
                w.writerow(row)
                f.flush()        # 每个突变算完就落盘，中断不丢
                print(f"[{n}/{len(tasks)}] {key}  "
                      f"ddG={row['ddG_mean(REU)']}±{row['ddG_std(REU)']} "
                      f"REU", flush=True)

        # 收尾：凡是还没写过的突变都补一行，包括一条轨迹都没成的。
        # 输出行数与 selected.csv 永远一致，缺谁按 nstruct_done 排序即可看出
        for key in sorted(set(by_key) - written):
            vals = got.get(key, [])
            row = summarize(by_key[key], vals, args.trials)
            w.writerow(row)
            f.flush()
            print(f"⚠️ {key} 只完成 {len(vals)}/{args.nstruct} 条轨迹  "
                  f"ddG={row['ddG_mean(REU)'] or '无'}", flush=True)

    print(f'\n总墙钟 {(time.time() - t0) / 3600:.2f} 小时')
    if fails:
        print('失败统计（现场目录已保留，可进去查空 db3）:')
        for key, (cnt, _, d) in sorted(fails.items()):
            print(f'  {key:<24} {cnt}/{args.nstruct} 条失败    {d}')
    print('输出:', args.out)


if __name__ == '__main__':
    main()
