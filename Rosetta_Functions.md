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

> **01 快速计算结合能 ΔG（不做 relax 与 repack）**

待补充 —— 脚本 `scripts/batch_dG_no_relax.py` 已完成，作为 02 的对照基线。直接对输入结构打分、沿 jump 拆开、直接相减，每个结构约 1 秒。

> **02 批量计算结合能 ΔG（relax + repack）**

把所有待评估的复合物 pdb 放进输入目录，批量扫完，每个结构输出一行指标。

核心是一个差值：**结合态的打分减去分离态的打分**，即 `dG_separated`。负得越多，说明结合让体系的能量降得越多。所谓「界面」在这里只是配角 —— 它决定 relax 的范围（省时间），以及 Rosetta 那个工具恰好叫 InterfaceAnalyzer。

参数集中在脚本顶部的配置块，日常只需改 `INTERFACE`：

```python
INPUTS    = '/data/lmk/rosetta_inputs'
OUTPUT    = '/data/lmk/rosetta_outputs/dG_relax_results.csv'
INTERFACE = 'HL_A'      # 链分组
RADIUS    = 8.0
```

```bash
python /data/lmk/rosetta_scripts/batch_dG_relax.py
```

每个结构的处理流程：

```
1. 读 PDB，自动补氢
2. 选出结合部位附近的残基                           决定 relax 范围，目前是CB < 8 Å
3. 只 relax 这些残基的侧链，主链固定
4. 重建 FoldTree，让 jump 1 恰好分开两组链
5. 打分
6. 沿 jump 把一组推开 → repack 侧链 → 打分
7. ΔG = E_complex − E_separated
```

**第 6 步的 repack 不能省**（`set_pack_separated(True)`）。两组分开后，原本朝向对方、被挤压的侧链会舒展开，这是真实的物理过程。保持结合态构象直接打分的话 `E_separated` 偏高，相减出来的 ΔG 被高估。


输出 csv 的列（`dG_relax_results.csv`）：

| 列 | 含义 |
| :--- | :--- |
| `pdb_id` | 文件名 |
| `nres` / `nres_relax` | 总残基数 / 参与 relax 的残基数 |
| `nres_antibody` / `nres_antigen` | 两组各自的残基数，**每次跑都该扫一眼**，对不上说明链分组有问题 |
| `E_complex(REU)` | 结合态的能量 |
| `E_antibody(REU)` / `E_antigen(REU)` | 分开并 repack 后，两部分各自的能量 |
| **`dG(REU)`** | **结合能，负得越多结合越强**；等于 `E_complex − E_antibody − E_antigen` |
| `dSASA_int(A^2)` | 界面埋藏面积 |
| `relax_time(s)` / `analyze_time(s)` | 两步各自耗时 |

⚠️ 列名把「下划线左边 = 抗体」写死了。`INTERFACE` 写成 `'A_HL'` 的话，`E_antibody` 里装的其实是抗原的能量，**而且不会报错**。

> **03 重复采样与统计**

待补充

> **04 多进程并行**

待补充

> **05 ΔΔG 突变扫描 / cartesian_ddG**

待补充

##### [PyRosetta API 文档](https://graylab.jhu.edu/PyRosetta.documentation/)
