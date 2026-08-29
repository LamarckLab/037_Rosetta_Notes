"""05 第一级：抗体侧饱和扫描。

跟 06 共用 ddG_screen_lib.py 里的同一套逻辑，只把默认值换成抗体侧，跑的时候不必
带参数。两个入口都做成薄封装而不是各存一份 —— 插入编号那个 bug 的教训是，重复的
代码意味着每个修复都要做两遍，而漏掉的那一遍不会有人发现。

判据：位点用重原子 4 Å（真实接触），repack 用邻居原子 8 Å（给侧链留响应空间）。
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import ddG_screen_lib as base

base.MUT_SIDE = 'antibody'
base.OUTPUT = '/data/lmk/rosetta_outputs/ddG_screen_antibody_results.csv'

if __name__ == '__main__':
    base.main()
