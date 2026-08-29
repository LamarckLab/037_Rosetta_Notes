"""06 第一级：抗原侧饱和扫描。

跟 05 共用 ddG_screen_lib.py 里的同一套逻辑，只把默认值换成抗原侧，跑的时候不必
带参数。做成薄封装而不是复制一份 —— 插入编号那个 bug 的教训是，重复的代码意味
着每个修复都要做两遍，而漏掉的那一遍不会有人发现。

判据与 05 相同：位点用重原子 4 Å（真实接触），repack 用邻居原子 8 Å（给侧链
留响应空间）。变的只有突变发生在哪一侧。
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import ddG_screen_lib as base

base.MUT_SIDE = 'antigen'
base.OUTPUT = '/data/lmk/rosetta_outputs/ddG_screen_antigen_results.csv'

if __name__ == '__main__':
    base.main()
