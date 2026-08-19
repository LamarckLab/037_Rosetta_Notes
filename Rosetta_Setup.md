## Lamarck &nbsp; &nbsp; &nbsp; 2026-08-05
#### 该文档用于部署 PyRosetta
---

## 01  创建 conda 环境
```bash
conda create -n lmk_Rosetta python=3.11 -y
conda activate lmk_Rosetta
```

## 02  安装依赖
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
> 脚本目录：/data/lmk/rosetta_tutorial  &nbsp;（.ipynb 与 .py）

## 06  在 conda 环境中配置并远程访问 JupyterLab
> 服务器与内核是解耦的：`jupyterlab` 提供服务器，`ipykernel` 让某个环境注册成为可选内核，一个 server 可挂多个内核，在网页上切换，新建一个环境 & 内核需要做

**(1) 安装 jupyterlab —— 使该环境能启动 jupyter 服务**
```bash
python -m pip install jupyterlab
```

**(2) 安装 ipykernel —— 使该环境能作为 jupyter 的 python 内核**
```bash
python -m pip install ipykernel
```

**(3) 把 lmk_Rosetta 环境注册为 jupyter 内核**
```bash
python -m ipykernel install --user --name lmk_rosetta --display-name lmk_Rosetta
```

**(4) 启动 jupyter 服务**
```bash
python -m jupyter lab --no-browser --ip=0.0.0.0 --port=8888
```

**(5) 在浏览器中连接 jupyter 前端**
```
http://192.168.208.236:8888/lab?token=<启动时打印的 token>
```
> 启动日志打印的是 `http://amax:8888/...`，该主机名在别的机器上解析不了，换成 IP 即可

##### [PyRosetta 官方下载页](https://www.pyrosetta.org/downloads)
