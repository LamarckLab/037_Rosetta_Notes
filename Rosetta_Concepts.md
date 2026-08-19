## Lamarck &nbsp; &nbsp; &nbsp; 2026-08-05
#### 该文档用于记录 PyRosetta 的核心概念，是看懂后续所有代码的前提
---


> **01 Pose**

Pose 是 Rosetta 表示「一个分子体系当前状态」的中心对象，所有操作都围绕它展开。里面装着构象（原子 xyz + 二面角 φ/ψ/χ）、序列与化学（每个位置是哪种残基、有哪些原子、如何连接）、链拓扑、上次打分的能量（含逐残基分解）、以及 PDBInfo（原始 PDB 的链号与残基编号）。

**Pose 是可变的**：需要留底用 `pose.clone()`。

不依赖文件也能建一个 Pose，用来单独观察这个对象：
```python
pose = pyrosetta.pose_from_sequence('AAAGGGKKK')
pose.total_residue()        # 9
pose.residue(3).name3()     # 'ALA'，残基编号从 1 开始，不是 0
pose.phi(3), pose.psi(3)    # (180.0, 180.0)，完全伸展
```
φ/ψ 全为 180° 说明 Pose 只是容器，构象是另外填进去的东西。

**两套残基编号** —— 最容易静默出错的地方。Rosetta 内部从 1 连续编到 N、无视链边界，PDB 文件则按链分别计数：

| 链   | Rosetta 编号 | PDB 编号 |
| :--- | :----------- | :------- |
| H    | 1 - 121      | 1 - 121  |
| L    | 122 - 230    | 1 - 109  |
| A    | 231 - 424    | 1 - 194  |

> **02 ScoreFunction**

打分函数：输入一个 Pose，输出一个数（REU）。构象好不好、突变有没有改善、binder 值不值得做，最终都归结为比较这个数。

```
总分 = Σ（权重ᵢ × 能量项ᵢ）
```

REF2015 共 19 项，半物理半统计 —— 一部分来自物理公式（`fa_atr` / `fa_rep` 范德华、`fa_elec` 静电、`fa_sol` 去溶剂化），一部分来自 PDB 数据库统计（`fa_dun` 侧链 rotamer 频率、`rama_prepro` 主链二面角分布）。

**REF2015** = Rosetta Energy Function 2015，当前默认的全原子能量函数，对应论文即本文末尾的 Alford et al. 2017。

**两个性质**

1. **单位是 REU**，与 kcal/mol 数量级接近但未经标定 —— 不要写成 kcal/mol，也不要拿去换算 Kd。
2. **是广延量**，蛋白越大总分越负。所以总分不能跨体系比较，评估 binder 必须用「分离前后的差值」而不是总分。

> **03 Mover**

待补充

> **04 ResidueSelector**

待补充

> **05 centroid 与 full-atom 两种精度表示**

待补充

##### [The Rosetta All-Atom Energy Function (Alford et al., 2017)](https://pubs.acs.org/doi/10.1021/acs.jctc.7b00125)
