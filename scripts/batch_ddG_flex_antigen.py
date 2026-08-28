"""06 第二级：对抗原侧突变跑官方 flex ddG。

跟 05 的 batch_ddG_flex.py 是同一套逻辑，只换默认的输入输出路径，跑的时候不必
带参数。协议本身完全相同 —— flex ddG 按 resfile 定位残基，不关心突变落在抗体
还是抗原，bubble 仍是官方原版的邻居原子 8 Å。

`--move-chain` 也不用改：它定义的是「算解离态时把哪组链移开」，与突变在哪一侧
无关，移开任一侧得到的都是同一个界面的解离。
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import batch_ddG_flex as base

base.SELECTED = '/data/lmk/rosetta_outputs/ddG_selected_antigen.csv'
base.OUTPUT = '/data/lmk/rosetta_outputs/ddG_flex_antigen_results.csv'

if __name__ == '__main__':
    base.main()
