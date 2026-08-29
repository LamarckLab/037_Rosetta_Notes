## Lamarck &nbsp; &nbsp; &nbsp; 2026-08-27
#### 该文档用于记录在 201server 上跑 PyRosetta 的命令
---

*201 机子的环境*
```bash
conda activate lmk_Rosetta
```

*输入输出路径*
```bash
输入目录:   /home/limingkai/rosetta_inputs
输出目录:   /home/limingkai/rosetta_outputs       # 指标 csv
脚本目录:   /home/limingkai/rosetta_scripts       # 生产脚本 .py
临时目录:   /home/limingkai/rosetta_work          # 只有 flex ddG 用，脚本自建
```

*并行线程*
```bash
--nproc 32      # 常规值；机子空闲时最多开到 64，留余量给同机其他人
```

⚠️ **本机 `nproc` 返回 176，但实际只有 82 个在线核心**。另外本机 **swap 为 0**，内存打满是直接 OOM 杀进程，没有变慢作为预警。

---

各档的原理、参数含义、输出 csv 的列，全部见 [Rosetta_Functions.md](./Rosetta_Functions.md)，本文档只放命令。

**命令一律显式带 `--inputs` / `--out`**：脚本配置块里写死的是 236 的 `/data/lmk/...`，用命令行覆盖掉，就不必在 201 上改脚本 —— 两台机子共用 git 里的同一份脚本，改动只发生在一处。

`--interface` 的左边必须是抗体、右边是抗原，且要与 pdb 实际链号一致，每批结构都先确认一遍。

所有脚本都**支持断点续跑**，中断后原样重跑即可。跳过的单位两档不同：01–04 按 `pdb_id`，05/06 按单个突变（`pdb_id` + 链 + 位置 + 插入编号 + 目标氨基酸）。

---

> **01 计算 ΔG：relax_no + repack_no**

```bash
python /home/limingkai/rosetta_scripts/batch_dG.py \
  --inputs /home/limingkai/rosetta_inputs \
  --out /home/limingkai/rosetta_outputs/dG_results.csv \
  --interface HL_A
```

> **02 计算 ΔG：relax_no + repack_yes**

```bash
python /home/limingkai/rosetta_scripts/batch_dG_repack.py \
  --inputs /home/limingkai/rosetta_inputs \
  --out /home/limingkai/rosetta_outputs/dG_repack_results.csv \
  --interface HL_A
```

> **03 计算 ΔG：relax_yes + repack_yes**

```bash
python /home/limingkai/rosetta_scripts/batch_dG_repack_relax.py \
  --inputs /home/limingkai/rosetta_inputs \
  --out /home/limingkai/rosetta_outputs/dG_repack_relax_results.csv \
  --interface HL_A
```

> **04 重复采样与并行**

```bash
python /home/limingkai/rosetta_scripts/batch_dG_pipeline.py \
  --inputs /home/limingkai/rosetta_inputs \
  --out /home/limingkai/rosetta_outputs/dG_pipeline_results.csv \
  --interface HL_A \
  --nstruct 5 \
  --nproc 32
```

> **05 ΔΔG 抗体侧饱和突变**

第一级。不用带 `--side`，入口本身就决定了突变哪一侧：

```bash
python /home/limingkai/rosetta_scripts/batch_ddG_screen_antibody.py \
  --inputs /home/limingkai/rosetta_inputs \
  --out /home/limingkai/rosetta_outputs/ddG_screen_antibody_results.csv \
  --interface HL_A \
  --nproc 32
```

跑完人工筛选，删掉不合适的行，另存为 `/home/limingkai/rosetta_outputs/ddG_selected_antibody.csv`，再进第二级：

```bash
python /home/limingkai/rosetta_scripts/batch_ddG_flex_antibody.py \
  --selected /home/limingkai/rosetta_outputs/ddG_selected_antibody.csv \
  --pdb-dir /home/limingkai/rosetta_inputs \
  --out /home/limingkai/rosetta_outputs/ddG_flex_antibody_results.csv \
  --xml /home/limingkai/rosetta_scripts/ddG-backrub.xml \
  --work /home/limingkai/rosetta_work/flexddg \
  --move-chain A \
  --nstruct 35 \
  --nproc 32
```

`--work` 目录脚本自己建，每条轨迹一个子目录，**成功失败都保留**（单条约 1.6 MB，一次跑几百条也就 1 GB 出头）。官方 nstruct=35 时单个突变约 18 CPU 小时，务必挂后台跑。

> **06 ΔΔG 抗原侧饱和突变**

跟 05 同一套逻辑与判据，只是突变发生在抗原侧，命令形状完全一样：

```bash
python /home/limingkai/rosetta_scripts/batch_ddG_screen_antigen.py \
  --inputs /home/limingkai/rosetta_inputs \
  --out /home/limingkai/rosetta_outputs/ddG_screen_antigen_results.csv \
  --interface HL_A \
  --nproc 32
```

人工筛选后另存为 `ddG_selected_antigen.csv`，再跑第二级：

```bash
python /home/limingkai/rosetta_scripts/batch_ddG_flex_antigen.py \
  --selected /home/limingkai/rosetta_outputs/ddG_selected_antigen.csv \
  --pdb-dir /home/limingkai/rosetta_inputs \
  --out /home/limingkai/rosetta_outputs/ddG_flex_antigen_results.csv \
  --xml /home/limingkai/rosetta_scripts/ddG-backrub.xml \
  --work /home/limingkai/rosetta_work/flexddg \
  --move-chain A \
  --nstruct 35 \
  --nproc 32
```

两侧的输出文件名都带侧别，**同时跑不会互相覆盖**；`--work` 可以共用，每条轨迹的子目录名里带了突变标识。

> **07 长任务挂后台**

201 是走公网的远程机子，SSH 一断前台任务就没了。04 与 05/06 的第二级动辄几小时到几天，一律挂后台：

```bash
nohup python /home/limingkai/rosetta_scripts/batch_ddG_flex_antibody.py \
  --selected /home/limingkai/rosetta_outputs/ddG_selected_antibody.csv \
  --pdb-dir /home/limingkai/rosetta_inputs \
  --out /home/limingkai/rosetta_outputs/ddG_flex_antibody_results.csv \
  --xml /home/limingkai/rosetta_scripts/ddG-backrub.xml \
  --work /home/limingkai/rosetta_work/flexddg \
  --nproc 32 \
  > /home/limingkai/rosetta_outputs/ddG_flex_antibody.log 2>&1 &
```

```bash
tail -f /home/limingkai/rosetta_outputs/ddG_flex_antibody.log      # 看进度
pgrep -af batch_ddG_flex                                  # 确认还在跑
```

##### [Rosetta_Functions.md](./Rosetta_Functions.md) &nbsp;|&nbsp; [Rosetta_Setup.md](./Rosetta_Setup.md)
