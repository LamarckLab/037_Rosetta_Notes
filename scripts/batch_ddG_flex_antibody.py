"""05 第二级：对抗体侧突变跑官方 flex ddG。

跟 06 共用 ddG_flex_lib.py 里的同一套逻辑，只换默认的输入输出路径，跑的时候不必
带参数。协议本身两侧完全相同 —— flex ddG 按 resfile 定位残基，不关心突变落在抗体
还是抗原，bubble 仍是官方原版的邻居原子 8 Å。
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import ddG_flex_lib as base

base.SELECTED = '/data/lmk/rosetta_outputs/ddG_selected_antibody.csv'
base.OUTPUT = '/data/lmk/rosetta_outputs/ddG_flex_antibody_results.csv'

if __name__ == '__main__':
    base.main()
