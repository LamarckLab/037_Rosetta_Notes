## Lamarck &nbsp; &nbsp; &nbsp; 2026-08-20
#### 该文档用于记录 server 上跑 PyRosetta 做批量评估的生产用法
---

*236 机子的环境*
```bash
conda activate lmk_Rosetta
```

*输入输出路径*
```bash
输入目录:   /data/lmk/rosetta_inputs       # 待评估的复合物 pdb，一个目录扫全部
输出目录:   /data/lmk/rosetta_outputs      # 指标 csv
脚本目录:   /data/lmk/rosetta_scripts      # 生产脚本 .py
```

---

**CSV 里有两列 ΔG，算的是同一个量，走的是两条独立的路**：

| 列            | 计算方式                                                         |
| :------------ | :--------------------------------------------------------------- |
| **`dG(REU)`** | 脚本自己算：把 Pose **真拆成两个**，可以拿到抗原和抗体各自的能量 |
| `dG_IA(REU)`  | InterfaceAnalyzerMover 算的：沿 jump 把一组链**平移**到远处      |

```
dG = E_complex − E_antibody − E_antigen
```

**三档工作流的分界线，即打分之前对侧链做了多少优化**：

|        | 侧链处理                | IA 的 `pack_input` / `pack_separated` | 耗时   |
| :----- | :---------------------- | :------------------------------------ | :----- |
| **01** | 什么都不做              | False / False                         | 秒级   |
| **02** | 界面 repack 一次        | True / True                           | 秒级   |
| **03** | 界面 FastRelax + repack | True / True                           | 分钟级 |

三档**都不动主链**，CA RMSD 恒为 0，结构不会偏离输入。**IA 本身不做 relax，只做 repack** —— 主链、键长、键角一概不动，只重新挑选侧链 rotamer。01 两条路都不做任何优化，所以两列 ΔG **必然完全相同**。

### 界面判据：两套，回答两个不同的问题

|                               | 判据   |    距离 | 测量对象        | 依据                             |
| :---------------------------- | :----- | ------: | :-------------- | :------------------------------- |
| **界面位点**                  | 重原子 | **4 Å** | 抗体组 ↔ 抗原组 | 结构生物学判定界面接触的通行口径 |
| **repack**                    | β 碳   | **8 Å** | 抗体组 ↔ 抗原组 | 沿用官方 flex ddG 的尺度         |
| **bubble**（仅 05/06 第二级） | β 碳   | **8 Å** | 突变位点 ↔ 周围 | 官方 flex ddG 原版，未改         |

**界面位点判据**回答「哪些残基真的在界面上」。01–04 用它填 `residues_num_interface` 这一列，05/06 用它决定枚举哪些位点做饱和突变。

**repack 判据**回答「算能量时放开哪些侧链去松弛」，判据要放宽，宁可多放几个也别漏掉能响应的侧链。

---

> **01 计算 ΔG：relax_no + repack_no**

把结构原样拿来打分，不做任何构象优化。**这不是领域标准做法**，官方对非 Rosetta 来源的结构建议先做侧链 repack。01 的用途是当诊断基线：跟 02 一比，就知道 repack 改变了多少；跟 03 一比，就知道 relax 又额外改变了多少。

```bash
python /data/lmk/rosetta_scripts/batch_dG.py
```

流程只有四步：

```
1. 读 PDB，自动补氢
2. 打分                                → E_complex
3. 按链拆成两个独立 Pose，各自打分       → E_antibody / E_antigen
4. ΔG = E_complex − E_antibody − E_antigen
```

**拆分用 `split_by_chain()` 而不是沿 jump 平移**，这样两部分的能量能分别拿到。

### 什么时候用

| 场景                 | 说明                                       |
| :------------------- | :----------------------------------------- |
| 快速过一遍几百个结构 | 每个约 2 秒                                |
| 判断结构本身干不干净 | `E_complex` 特别大说明输入结构带着严重应变 |
| 作为 02 / 03 的对照  | 差值分别是 repack 与 relax 的贡献          |

> **02 计算 ΔG：relax_no + repack_yes**

**Rosetta 官方推荐的标准做法。** 不做 relax，但结合态与分离态各做一次 repack。官方文档对非 Rosetta 来源的结构（AF3 / RFdiffusion 输出等）正是这么建议的：

```bash
python /data/lmk/rosetta_scripts/batch_dG_repack.py
```

流程五步：

```
1. 读 PDB，自动补氢
2. 选出结合部位附近的残基                     决定 repack 范围，CB < 8 Å
3. 结合态 repack 这批残基，打分                → E_complex
4. 按链拆成两个 Pose，各自 repack 同一批残基    → E_antibody / E_antigen
5. ΔG = E_complex − E_antibody − E_antigen
```

### 什么时候用

| 场景                   | 说明                                     |
| :--------------------- | :--------------------------------------- |
| **默认选它**           | 官方标准，比 01 可靠、比 03 快 10 倍左右 |
| 结构来自晶体或预测模型 | 正是官方建议先做一次 repack 的场景       |

02 只做**一次离散 repack**，侧链只能落在 rotamer 库里现成的位置上。输入结构应变严重时（`E_complex` 极大）这一步往往不够，需要 03 的多轮迭代加连续最小化。

### `dG` 与 `dG_IA` 为什么在这一档开始分歧

01 里两条路都不做任何优化，两列必然相同。打开 repack 之后才出现差异：

| 步骤            | `dG`（脚本自算）                       | `dG_IA`（InterfaceAnalyzer）         |
| :-------------- | :------------------------------------- | :----------------------------------- |
| 起点            | 原始 pose                              | 同一个 pose 的 clone（**起点相同**） |
| repack 哪些残基 | 界面位点 CB 8 Å 选出的那批             | 界面位点 CB 8 Å 选出的同一批         |
| 结合态          | repack 那批 → 打分                     | `pack_input=True` 自己 repack → 打分 |
| **怎么分开**    | **`split_by_chain()` 真拆成两个 Pose** | **沿 jump 平移到远处**               |
| 分离态          | 两个子 Pose 各 repack 同一批 → 打分    | `pack_separated=True` 再 repack      |

界面判据对齐之后，两条路的差异只剩**分开方式**：脚本真拆成两个 Pose，IA 沿 jump 平移到远处。

> **03 计算 ΔG：relax_yes + repack_yes**

把 02 的「一次 repack」换成 FastRelax。**自由度跟 02 完全相同** —— 主链全程固定，只动界面那批侧链 —— 区别在于搜索方式：FastRelax 是多轮 repack + 连续最小化，过程中给 `fa_rep` 加斜坡，让纠缠的原子先错开再收紧。彻底得多，代价是慢十倍。

```bash
python /data/lmk/rosetta_scripts/batch_dG_repack_relax.py
```

流程七步：

```
1. 读 PDB，自动补氢
2. 选出结合部位附近的残基                     决定 relax 与 repack 范围，CB < 8 Å
3. 只 relax 这批残基的侧链，主链固定
4. 结合态 repack 同一批残基，打分               → E_complex
5. 按链拆成两个 Pose，重建二硫键
6. 各自 repack 同一批残基，再分别打分           → E_antibody / E_antigen
7. ΔG = E_complex − E_antibody − E_antigen
```

### 什么时候用

| 场景             | 说明                                               |
| :--------------- | :------------------------------------------------- |
| 输入结构应变严重 | `E_complex` 极大时，一次 repack 跳不出坏的局部极小 |
| 需要更稳健的结果 | 连续最小化能跳出 rotamer 库的离散限制              |

### 输出 csv 的列

**02 与 03 的列完全一致**，01 少一列 `residues_num_interface`（它不做残基选择）。

| 列                                               | 含义                                                                    |
| :----------------------------------------------- | :---------------------------------------------------------------------- |
| `pdb_id`                                         | 文件名                                                                  |
| `residues_num_total`                             | 复合物的总残基数                                                        |
| `residues_num_interface`                         | 真实界面大小，重原子 4 Å 判据                                           |
| `residues_num_repack`                            | 实际 repack 的残基数，CB 8 Å 判据（01 无此列，它不 repack）             |
| `residues_num_antibody` / `residues_num_antigen` | 抗体、抗原各自的残基数                                                  |
| `E_complex(REU)`                                 | 结合态的能量                                                            |
| `E_antibody(REU)` / `E_antigen(REU)`             | 分开后两部分各自的能量（02 / 03 会先 repack）                           |
| **`dG(REU)`**                                    | **结合能，负得越多结合越强**；等于 `E_complex − E_antibody − E_antigen` |
| `dG_IA(REU)`                                     | InterfaceAnalyzer 独立算的 ΔG，供对照；两列越接近，该结构的 ΔG 越可信   |
| `dSASA_int(A^2)`                                 | 分开前后溶剂可及面积之差，即界面埋藏的面积                              |
| `total_time(s)`                                  | 该结构从读 PDB 到算完的总耗时                                           |

> **04 重复采样与并行**

生产级别脚本：batch_dG_pipeline.py

在 03 的基础上做两件事：**每个结构重复 `NSTRUCT` 遍取统计量**，以及**按结构切分到 `NPROC` 个进程上并行**。

```bash
python /data/lmk/rosetta_scripts/batch_dG_pipeline.py
```

**`dG` 各列取的是 InterfaceAnalyzer 的值**，即 03 CSV 里的 `dG_IA` 那条路，不是 03 的 `dG` 列。IA 是业界通用实现，所以这里不再加 `_IA` 后缀。之前两条路都留，是为了拿到 `E_antibody` / `E_antigen`。

PyRosetta 每个进程自动取随机种子（`run:constant_seed = False`），同进程内连续调用也在推进 RNG，所以不需要手动管种子。

### `--nproc` 怎么定

Rosetta 是**单线程 CPU 密集型**，一个进程占一个核心，进程数直接等于占用的核心数。实际取 `min(--nproc, 待算结构数)`。

### 输出 csv 的列

| 列                            | 含义                                     |
| :---------------------------- | :--------------------------------------- |
| `pdb_id`                      | 文件名                                   |
| `residues_num_*`              | 与 03 相同的四列残基数                   |
| `nstruct`                     | 实际重复了几遍                           |
| **`dG_mean(REU)`**            | **ΔG 均值，排序看这一列**（来自 IA）     |
| **`dG_std(REU)`**             | **标准差，两个设计的差距要大于它才可信** |
| `dG_min(REU)` / `dG_max(REU)` | 采样区间，一眼看出最坏情况               |
| `dG_all(REU)`                 | 每一遍的原始值，分号分隔，留着重新统计   |
| `dSASA_int_mean(A^2)`         | 界面埋藏面积的均值                       |
| `total_time(s)`               | 该结构 `nstruct` 遍加起来的总耗时        |

> **05 ΔΔG：抗体侧饱和突变**

ΔG 回答「这个复合物结合得多牢」，ΔΔG 回答「**改一个氨基酸会让它更牢还是更松**」——亲和力成熟要的是后者。

```
ΔΔG = ΔG_突变体 − ΔG_野生型          ΔΔG < 0 表示突变体结合更强
```

### 为什么分两级

界面结合 ΔΔG 的领域标准是 **flex ddG**，用 backrub 让主链微动。它准，但贵：

```
单条轨迹 35000 backrub steps   ≈ 31 分钟
官方 nstruct = 35              ≈ 18 CPU 小时 / 单个突变
```

一个界面饱和扫描是几百个突变，用 flex ddG 就是 **32 核跑一个月**。所以拆成两级：

|            | 用什么                       | 速度   | 产出             |
| :--------- | :--------------------------- | :----- | :--------------- |
| **第一级** | REF2015 + 只 repack（即 02） | 秒级   | **只是排序**     |
| **第二级** | 官方 flex ddG（talaris2014） | 小时级 | **可报告的 ΔΔG** |

中间是**人工筛选**：第一级把全部突变排好序输出，删掉不合适的行另存一份，第二级只算留下的。

---

#### 第一级：饱和扫描

```bash
python /data/lmk/rosetta_scripts/batch_ddG_screen_antibody.py
```

不用准备突变列表。脚本取界面**抗体侧**的残基（重原子 4 Å 判据），每个位点枚举其余 19 种氨基酸。

```
1. 选界面残基 —— 两侧全要，CB 8 Å           → repack 范围
2. 选可突变位点 —— 只要一侧，重原子 4 Å      → 突变范围
3. 算野生型                                  → dG_wt
4. 逐个突变：MutateResidue → 按 02 的流程算  → dG_mut
5. ddG_screen = dG_mut − dG_wt
6. 全部算完后按 ddG_screen 升序重排，补 rank
```

**突变范围 ≠ repack 范围**：只突变抗体侧，但两侧都 repack。抗原侧虽不突变，但会 repack。

**`dG_wt` 要抽样多次再截尾平均**（脚本定义的是 20 次去掉两端各 5 个）。

### ⚠️ 插入编号：抗体必踩的坑

Kabat / Chothia 编号里 CDR 常有 **H100A、H52A** 这种插入编号 —— 它和 H100、H52 是**不同的残基**。

所以 csv 里 `pdb_position` 之外还有一列 **`icode`**，两列合起来才唯一定位一个残基。

### ⚠️ 这一列为什么叫 `ddG_screen` 而不是 `ddG`

只 repack 没有主链让步的余地，导致已实测的系统性偏差：

```
S → Y 给出 1405 REU 这种物理上无意义的数
```

所以它**只能用于排序，不能作为 ΔΔG 报告**。报告值一律来自第二级。

---

#### 第二级：官方 flex ddG

```bash
python /data/lmk/rosetta_scripts/batch_ddG_flex_antibody.py
```

筛选后的 csv（`ddG_selected_antibody.csv`），跑 Kortemme 实验室的 `ddG-backrub.xml`。该 XML 原样存在 `scripts/` 下。

```
1. 按 位置/链/目标氨基酸 生成 resfile
2. 跑 NSTRUCT 条独立的 backrub 轨迹
3. 每条轨迹取四个状态的逐项能量
4. ΔΔG = (bound_mut + unbound_wt) − (bound_wt + unbound_mut)
5. 过 ZEMu GAM 重加权 → kcal/mol
6. 对 NSTRUCT 条求均值与标准差
```

**任务粒度是一条轨迹，不是一个突变。** 按突变切的话每个任务几十小时，进程之间没法均衡。

### GAM 重加权

第 4 步相减得到的是各能量项的 ΔΔG，单位是 REU，直接加起来与实验值对不上。flex ddG 的做法是把 7 个分子间项（`fa_atr`、`fa_rep`、`fa_sol`、`fa_elec` 与三个氢键项）各自过一条拟合好的 S 形曲线再求和 —— 这就是 **GAM**（generalized additive model，广义加性模型），系数由作者在 ZEMu 数据集上拟合。

```
ddG_mean(kcal/mol)     过了 GAM，标定到 kcal/mol，与实验可比  ← 报告用这个
ddG_raw_mean(REU)      7 项直接相加，未标定                    ← 只作对照
```

论文报的 MAE ≈ 1 kcal/mol 是 GAM 之后的值。两列的符号通常一致，差得远说明这个突变落在拟合曲线的非线性区，结果要打折扣看。

### ⚠️ 本档用 talaris2014

flex ddG 整套是在 **talaris2014** 上标定的，与 01–04 不可比。

### 输出 csv 的列

| 列                                            | 含义                                            |
| :-------------------------------------------- | :---------------------------------------------- |
| `pdb_id` / `chain` / `pdb_position` / `icode` | 突变位置；`icode` 是插入编号，多数为空          |
| `wt_aa` / `mut_aa`                            | 原氨基酸 / 目标氨基酸                           |
| **`ddG_mean(kcal/mol)`**                      | **业界报告值**，GAM 重加权后对 nstruct 求均值   |
| `ddG_std(kcal/mol)`                           | nstruct 条轨迹之间的标准差                      |
| `ddG_all(kcal/mol)`                           | 每条轨迹的原始值，分号分隔，换统计口径不必重跑  |
| `ddG_raw_mean(REU)` / `ddG_raw_std(REU)`      | 未经 GAM 重加权的 talaris2014 加和，供对照      |
| `nstruct_done`                                | 实际完成几条轨迹；为 0 表示该突变一条也没取到数 |
| `backrub_trials`                              | 每条轨迹的 backrub 步数                         |
| `cpu_time(s)`                                 | 各条轨迹耗时之和，即该突变真实占用的 CPU 时间   |

> **06 ΔΔG：抗原侧饱和突变**

跟 05 **完全同一套逻辑与判据**，只是突变发生在抗原侧。脚本是独立入口，跑的时候不用带参数：

```bash
# 第一级
python /data/lmk/rosetta_scripts/batch_ddG_screen_antigen.py

# 人工筛选，另存为 ddG_selected_antigen.csv

# 第二级
python /data/lmk/rosetta_scripts/batch_ddG_flex_antigen.py
```

### 与 05 的差别只有一处

|               | 05                                                 | 06                                                                 |
| :------------ | :------------------------------------------------- | :----------------------------------------------------------------- |
| 突变哪一侧    | 抗体（`HL_A` 下划线左边）                          | 抗原（`HL_A` 下划线右边）                                          |
| 位点判据      | 重原子 4 Å                                         | 同                                                                 |
| repack 判据   | CB 8 Å                                             | 同                                                                 |
| 第二级 bubble | 官方 CB 8 Å                                        | 同                                                                 |
| 输出文件      | `ddG_screen_antibody_results.csv`<br>`ddG_flex_antibody_results.csv` | `ddG_screen_antigen_results.csv`<br>`ddG_flex_antigen_results.csv` |

第二级的**协议本身完全没变** —— flex ddG 按 resfile 定位残基，不关心突变落在抗体还是抗原。`--move-chain` 也不用改：它定义的是「算解离态时移开哪组链」，与突变在哪一侧无关，移开任一侧得到的都是同一个界面的解离。

抗原侧的位点通常比抗体侧少：实测同一结构抗体侧 17 个、抗原侧 15 个。

---

##### [flex ddG 教程与协议](https://github.com/Kortemme-Lab/flex_ddG_tutorial) &nbsp;|&nbsp; [Barlow et al. JPCB 2018](https://pubs.acs.org/doi/10.1021/acs.jpcb.7b11367)

##### [PyRosetta API 文档](https://graylab.jhu.edu/PyRosetta.documentation/)
