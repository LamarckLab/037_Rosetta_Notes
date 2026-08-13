<h1 align="center">🧬 Rosetta 部署与功能测试笔记</h1>

<p align="center"><em>—— 2026.08.05</em></p>

<p align="center">
  <img src="https://img.shields.io/badge/Tool-Rosetta-blue?style=flat-square" />
  <img src="https://img.shields.io/badge/Field-Protein%20Design-orange?style=flat-square" />
  <img src="https://img.shields.io/badge/Platform-Linux-555?style=flat-square" />
  <img src="https://img.shields.io/badge/Status-In%20Progress-yellow?style=flat-square" />
</p>

---

## 内容索引

| 文档                                                     | 说明                                                                     |
| :------------------------------------------------------- | :------------------------------------------------------------------------ |
| [Rosetta_Setup.md](./Rosetta_Setup.md)                   | conda 环境、依赖与 PyRosetta 安装、验证、工作目录、JupyterLab 启动        |
| [Rosetta_Concepts.md](./Rosetta_Concepts.md)             | 五个核心概念：Pose / ScoreFunction / Mover / ResidueSelector / 精度表示   |
| [Rosetta_Functions.md](./Rosetta_Functions.md)           | 各功能代码：结构打分 / FastRelax / 界面分析 / 逐残基能量 / ddG 扫描 / 批量并行 |
| [Rosetta_Metrics_Format.md](./Rosetta_Metrics_Format.md) | InterfaceAnalyzer 输出字段参考                                           |
| [input_data/](./input_data/)                             | 各功能测试用输入文件                                                     |

---

##### [Rosetta 官方文档](https://docs.rosettacommons.org) &nbsp;|&nbsp; [PyRosetta 官网](https://www.pyrosetta.org)
