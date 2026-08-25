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

| 列            | 怎么来的                                                         |
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

**界面判据全仓库统一为：跨链两个残基的 CB 距离 ≤ 10 Å**（甘氨酸没有 CB，用 CA）。这是 InterfaceAnalyzer 的默认值。

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
2. 选出结合部位附近的残基                     决定 repack 范围，CB < 10 Å
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
| repack 哪些残基 | CB 10 Å 选出的那批                     | CB 10 Å 选出的同一批                 |
| 结合态          | repack 那批 → 打分                     | `pack_input=True` 自己 repack → 打分 |
| **怎么分开**    | **`split_by_chain()` 真拆成两个 Pose** | **沿 jump 平移到远处**               |
| 分离态          | 两个子 Pose 各 repack 同一批 → 打分    | `pack_separated=True` 再 repack      |

界面判据对齐之后，主要只剩两条差异：**分开方式**（真拆成两个 Pose vs 平移到远处），以及 **repack 本身是模拟退火**，同一批残基跑两遍也不会落在完全相同的解上。

所以 `dG_IA` 这一列的作用是**可靠性探针**：两个数越接近，说明该结构的 ΔG 越不受实现细节影响。分歧大到会改变设计之间的排序时，用 04 的重复采样取均值。

> **03 计算 ΔG：relax_yes + repack_yes**

把 02 的「一次 repack」换成 FastRelax。**自由度跟 02 完全相同** —— 主链全程固定，只动界面那批侧链 —— 区别在于搜索方式：FastRelax 是多轮 repack + 连续最小化，过程中给 `fa_rep` 加斜坡，让纠缠的原子先错开再收紧。彻底得多，代价是慢十倍。

```bash
python /data/lmk/rosetta_scripts/batch_dG_repack_relax.py
```

流程七步：

```
1. 读 PDB，自动补氢
2. 选出结合部位附近的残基                     决定 relax 与 repack 范围，CB < 10 Å
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
| `residues_num_interface`                         | 界面残基数，即 CB 10 Å 选出的那批（01 无此列）                          |
| `residues_num_antibody` / `residues_num_antigen` | 抗体、抗原各自的残基数                                                  |
| `E_complex(REU)`                                 | 结合态的能量                                                            |
| `E_antibody(REU)` / `E_antigen(REU)`             | 分开后两部分各自的能量（02 / 03 会先 repack）                           |
| **`dG(REU)`**                                    | **结合能，负得越多结合越强**；等于 `E_complex − E_antibody − E_antigen` |
| `dG_IA(REU)`                                     | InterfaceAnalyzer 独立算的 ΔG，供对照；两列越接近，该结构的 ΔG 越可信   |
| `dSASA_int(A^2)`                                 | 分开前后溶剂可及面积之差，即界面埋藏的面积                              |
| `total_time(s)`                                  | 该结构从读 PDB 到算完的总耗时                                           |

> **04 重复采样与统计**

待补充

> **05 多进程并行**

待补充

> **06 ΔΔG 突变扫描 / cartesian_ddG**

待补充

##### [PyRosetta API 文档](https://graylab.jhu.edu/PyRosetta.documentation/)
