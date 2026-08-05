## Lamarck &nbsp; &nbsp; &nbsp; 2026-08-05
#### 该文档用于部署 PyRosetta
---

## 01  创建 conda 环境
```bash
conda create -n lmk_Rosetta python=3.11 -y
conda activate lmk_Rosetta
```

## 02  安装依赖
> pandas 用于出指标表，ipykernel 供 Remote-SSH 里的 notebook 选到该环境
```bash
python -m pip install pandas ipykernel
```

## 03  安装 PyRosetta
> 非商业用途免费，无需申请 license；官方推荐直接 pip 装季度稳定版
```bash
python -m pip install pyrosetta --find-links https://west.rosettacommons.org/pyrosetta/quarterly/release
```

## 04  验证
> `-mute all` 静默 init 的上百行日志；
```bash
python -c "import pyrosetta; pyrosetta.init('-mute all'); print('PyRosetta OK')"
```

## 05  PyRosetta 输入输出目录
**236 机子路径**
> 输入目录：/data/lmk/rosetta_inputs  &nbsp;（待评估的复合物 pdb）
> 输出目录：/data/lmk/rosetta_outputs  &nbsp;（relax 后结构 + 指标 csv）

##### [PyRosetta 官方下载页](https://www.pyrosetta.org/downloads)
