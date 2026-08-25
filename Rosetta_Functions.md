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

**主指标是 `dG_IA`** —— InterfaceAnalyzerMover 算出来的结合能。这是权威的 dG 计算方法，不依赖脚本自己选的界面判据。CSV 里另一列 `dG` 是按照脚本定义的逻辑算的。

InterfaceAnalyzer 拿到 Pose 之后做四件事，其中两步由开关控制：

```
1. 若 pack_input=True        对复合物 repack        → E_complex
2. 沿 jump 把一组链平移到远处
3. 若 pack_separated=True    对分离态 repack        → E_separated
4. ΔG = E_complex − E_separated
```

**两个开关都是可配置的，这是三档工作流的分界线**：

|        | `pack_input` | `pack_separated` | 脚本另外做了什么         |
| :----- | :----------- | :--------------- | :----------------------- |
| **01** | False        | False            | 无                       |
| **02** | True         | True             | 无                       |
| **03** | True         | True             | 交给 IA 之前先 FastRelax |

⚠️ **IA 本身不做 relax，只做 repack** —— 主链、键长、键角一概不动，只重新挑选侧链 rotamer。

**`dG` 与 `dG_IA` 的差 = repack 环节的方法学差异。** 两者算的是同一个量，区别在于界面残基集合的定义方式：（脚本里用 CB 8 Å，IA 用分离前后 SASA 变化）、分开方式（脚本真拆成两个 Pose，IA 沿 jump 平移）。01 里两个 repack 开关都关闭，所以两列**必然完全相同**。02 和 03 里两者会有差异，差得大说明该结构的 ΔG 对界面定义敏感。

---

> **01 计算 ΔG：relax_no + repack_no**

把结构原样拿来打分，不做任何构象优化。**这不是领域标准做法**，官方对非 Rosetta 来源的结构建议开启两个 pack 开关。01 的用途是当诊断基线：跟 02 一比，就知道 repack 改变了多少；跟 03 一比，就知道 relax 又额外改变了多少。

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
| 结构来自晶体或预测模型 | 正是官方建议开两个 pack 开关的场景       |
| 批量筛几百上千个       | 每个约 5 秒                              |

02 只做**一次离散 repack**，侧链只能落在 rotamer 库里现成的位置上。输入结构应变严重时（`E_complex` 极大）这一步往往不够，需要 03 的多轮迭代加连续最小化。

### `dG` 与 `dG_IA` 为什么在这一档开始分歧

打开 repack 开关之后，两条路的做法出现四处差异：

| 步骤                | `dG`（脚本自算）                    | `dG_IA`（InterfaceAnalyzer）         |
| :------------------ | :---------------------------------- | :----------------------------------- |
| 起点                | 原始 pose                           | 同一个 pose 的 clone（**起点相同**） |
| **repack 哪些残基** | **CB 8 Å 选出的那批**               | **它自己判定的一批，范围明显更大**   |
| 结合态              | repack 那批 → 打分                  | `pack_input=True` 自己 repack → 打分 |
| 怎么分开            | `split_by_chain()` 真拆成两个 Pose  | 沿 jump 平移到远处                   |
| 分离态              | 两个子 Pose 各 repack 同一批 → 打分 | `pack_separated=True` 再 repack      |

**主因是 repack 范围。** 证据在耗时两列：`repack_time` 不到半秒，而 `analyze_time` 要 4~5 秒 —— 同样是 repack 两次，IA 慢一个数量级，说明它处理的残基多得多。
范围不同 → 消除的应变量不同 → 两个状态各自偏移不同的量 → 相减出来的 ΔG 就不一样。

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

三档**都是主链固定**，CA RMSD 恒为 0，结构不会偏离输入。

### 什么时候用

| 场景             | 说明                                               |
| :--------------- | :------------------------------------------------- |
| 输入结构应变严重 | `E_complex` 极大时，一次 repack 跳不出坏的局部极小 |
| 需要更稳健的结果 | relax 清掉应变后，`dG` 与 `dG_IA` 的分歧显著收窄   |
| 结构数量不多     | 每个约 1 分钟                                      |

### 输出 csv 的列

三档的列基本一致，01 没有 `nres_relax` / `relax_time(s)`，02 对应位置是 `nres_repack` / `repack_time(s)`。

| 列                                   | 含义                                                                         |
| :----------------------------------- | :--------------------------------------------------------------------------- |
| `pdb_id`                             | 文件名                                                                       |
| `nres` / `nres_relax`                | 总残基数 / 参与 relax 的残基数                                               |
| `nres_antibody` / `nres_antigen`     | 抗体抗原各自的残基数                                                         |
| `E_complex(REU)`                     | 结合态的能量                                                                 |
| `E_antibody(REU)` / `E_antigen(REU)` | 分开并 repack 后，两部分各自的能量                                           |
| **`dG(REU)`**                        | **结合能，负得越多结合越强**；等于 `E_complex − E_antibody − E_antigen`      |
| `dG_IA(REU)`                         | InterfaceAnalyzer 自己算的 ΔG，供对照；与 `dG` 的差 = 两边 repack 实现的差异 |
| `dSASA_int(A^2)`                     | 分开前后溶剂可及面积之差，即界面埋藏的面积                                                           |
| `relax_time(s)` / `analyze_time(s)`  | 两步各自耗时                                                                 |

> **04 重复采样与统计**

待补充

> **05 多进程并行**

待补充

> **06 ΔΔG 突变扫描 / cartesian_ddG**

待补充

##### [PyRosetta API 文档](https://graylab.jhu.edu/PyRosetta.documentation/)
