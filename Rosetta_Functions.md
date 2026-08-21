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

> **00 批量计算界面结合能**

把所有待评估的复合物 pdb 放进输入目录，一条命令扫完，每个结构输出一行指标。

参数集中在脚本顶部的配置块，日常只需改 `INTERFACE`：

```python
INPUTS    = '/data/lmk/rosetta_inputs'
OUTPUT    = '/data/lmk/rosetta_outputs/interface_metrics.csv'
INTERFACE = 'HL_A'      # 链分组，要与 pdb 里的实际链号一致
RADIUS    = 8.0
```

长任务放 screen 里跑，断线不中断：

```bash
screen -S batch
python /data/lmk/rosetta_scripts/batch_interface.py
```

每个结构的处理流程：

```
1. 读 PDB，自动补氢                                      PDB 通常只有重原子
2. 选界面残基                                            距对方链 CB < 8 Å 的两侧残基
3. 只 relax 这些残基的侧链，主链固定                       ← 约 50 s，耗时全在这
4. 重建 FoldTree，让 jump 1 恰好分开两组链
5. 打分                                                  → E_complex
6. 沿 jump 把一组推开 → repack 侧链 → 打分                → E_separated
7. ΔG = E_complex − E_separated                          ← 4~7 一次 apply 完成，约 3 s
```

**第 6 步的 repack 不能省**（`set_pack_separated(True)`）。两组分开后，原本朝向对方、被挤压的侧链会舒展开，这是真实的物理过程。保持结合态构象直接打分的话 `E_separated` 偏高，相减出来的 ΔG 被高估。

⚠️ **repack ≠ relax**：repack 只重新挑选侧链 rotamer，不动主链、不做梯度最小化，所以第 6 步只要几秒；第 3 步的 relax 是 repack + 最小化多轮迭代，要 50 秒。

**第 5~6 步不是真把结构切成两个文件**，而是沿 jump 把一组链平移到很远处，仍在同一个 Pose 里打分 —— 距离足够远，相互作用项自然归零。这也是第 4 步必须先重建 FoldTree 的原因：得先有一个恰好分开两组的 jump，才有东西可推。

**为什么只 relax 界面**：主链固定保证结构不偏离输入（CA RMSD = 0），只处理界面能把耗时压到全侧链 relax 的一半。剩余部分未处理带来的应变，在「复合物 − 分开的两部分」这个差值里会抵消。

**界面分组不能写死 jump 编号** —— 不同 pdb 的链顺序不同。脚本用 `setup_foldtree(pose, 'HL_A', jumps)` 按链分组重建 FoldTree，让 jump 1 恰好分开两组，再交给 `InterfaceAnalyzerMover(1)`：

```python
jumps = vector1_int()
setup_foldtree(pose, spec, jumps)     # 重建后 jump 1 = 两组之间
ia = InterfaceAnalyzerMover(1)
ia.set_pack_separated(True)           # 分开后重新 repack，模拟解离时侧链舒展
ia.set_calc_dSASA(True)               # dSASA 默认不算，要显式打开
```

⚠️ **链分组只能用字符串传给 `setup_foldtree`**。`InterfaceAnalyzerMover('HL_A')` 和 `ia.set_interface('HL_A')` 都不接受裸字符串，它们要的是 `DockingPartners` 对象，而该对象没有字符串构造入口。

输出 csv 的列：

| 列                    | 含义                                              |
| :-------------------- | :------------------------------------------------ |
| `pdb`                 | 文件名                                            |
| `nres` / `nres_iface` | 总残基数 / 界面残基数                             |
| **`dG_separated`**    | **界面结合能 REU，负得越多结合越强**              |
| `dSASA`               | 界面埋藏面积 Å²                                   |
| `total_score`         | relax 后的总分，仅供参考（局部 relax 后必然偏高） |
| `relax_s` / `ia_s`    | 两步各自耗时                                      |

> 断点续跑：每算完一个立即追加并 flush，重跑时自动跳过 csv 里已有的文件名。中断后直接再跑一次即可。
> 单个结构失败只打印错误继续下一个，不中断整批。

> **01 重复采样与统计**

待补充

> **02 多进程并行**

待补充

> **03 点突变扫描 / cartesian_ddG**

待补充

##### [PyRosetta API 文档](https://graylab.jhu.edu/PyRosetta.documentation/)
