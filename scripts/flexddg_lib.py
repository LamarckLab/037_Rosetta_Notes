"""flex ddG 的取数与重加权，照 Kortemme 实验室官方 analyze_flex_ddG.py 实现。

单独成一个模块，是因为这部分是纯数据处理、不依赖 PyRosetta，可以脱离服务器单测。

参考：
  Barlow et al. J. Phys. Chem. B 2018, "Flex ddG: Rosetta Ensemble-Based
  Estimation of Changes in Protein-Protein Binding Affinity upon Mutation"
  https://github.com/Kortemme-Lab/flex_ddG_tutorial
"""

import math
import sqlite3

# ZEMu GAM 的拟合参数，逐能量项一组 (a, b)。数值取自官方脚本，不要改
ZEMU_GAM_PARAMS = {
    'fa_sol':      (6.940, -6.722),
    'hbond_sc':    (1.902, -1.999),
    'hbond_bb_sc': (0.063,  0.452),
    'fa_rep':      (1.659, -0.836),
    'fa_elec':     (0.697, -0.122),
    'hbond_lr_bb': (2.738, -1.179),
    'fa_atr':      (2.313, -1.649),
}

# db3 里除各能量分量外还存了 total_score，它本身就是各分量之和，求和时必须排除
DERIVED_TERMS = {'total_score'}

# 四个状态：ΔΔG = (bound_mut + unbound_wt) - (bound_wt + unbound_mut)
STATE_SIGN = {'bound_mut': +1, 'unbound_wt': +1,
              'bound_wt': -1, 'unbound_mut': -1}


def gam(x, term):
    """把一个原始 REU 分量映射到 GAM 拟合后的尺度（一条缩放过的 sigmoid）"""
    a, b = ZEMU_GAM_PARAMS[term]
    return -math.exp(a) + 2.0 * math.exp(a) / (1.0 + math.exp(-x * math.exp(b)))


def ddG_from_terms(terms):
    """逐项 ΔΔG -> (GAM 重加权总和, 原始加和)

    terms: {score_term: 该项的 ΔΔG}，至少要含 ZEMU_GAM_PARAMS 的 7 项
    """
    missing = set(ZEMU_GAM_PARAMS) - set(terms)
    if missing:
        raise KeyError(f'db3 里缺少能量项: {sorted(missing)}')
    gam_total = sum(gam(terms[t], t) for t in ZEMU_GAM_PARAMS)
    raw_total = sum(v for k, v in terms.items() if k not in DERIVED_TERMS)
    return gam_total, raw_total


# InterfaceDdGMover 内部按固定顺序把四个状态交给 db_reporter，batch_id 依次为 1..4。
# 这个顺序在 C++ 里写死，XML 看不到；下面的 read_db3 会在 batches.name 可用时反查校验。
BATCH_ORDER = ('bound_wt', 'unbound_wt', 'bound_mut', 'unbound_mut')

# batches 用 LEFT JOIN：它偶尔是空表（结构在、无数据），而能量本身存在
# structure_scores 里并不受影响。用 INNER JOIN 会让整个查询返回零行。
_SCORE_SQL = """
SELECT structure_scores.batch_id, batches.name, structure_scores.struct_id,
       score_types.score_type_name, structure_scores.score_value
FROM structure_scores
INNER JOIN score_types
        ON structure_scores.score_type_id = score_types.score_type_id
       AND structure_scores.batch_id = score_types.batch_id
INNER JOIN structures ON structures.struct_id = structure_scores.struct_id
LEFT  JOIN batches    ON batches.batch_id = structures.batch_id
"""


def state_of(batch_id, batch_name):
    """定状态：优先用 batches.name，空表时退回 batch_id 顺序"""
    if batch_name:
        named = match_state(batch_name)
        by_id = BATCH_ORDER[batch_id - 1] if 1 <= batch_id <= 4 else None
        if named and by_id and named != by_id:
            raise ValueError(f'batch_id {batch_id} 的名字是 {batch_name}，'
                             f'与假定顺序 {by_id} 不符，BATCH_ORDER 需要更新')
        return named
    return BATCH_ORDER[batch_id - 1] if 1 <= batch_id <= 4 else None


def read_db3(path):
    """从一次 flex ddG 的 ddG.db3 里读出 {状态: {能量项: 值}}

    trajectory_stride 小于总 trials 时，一个状态会存多个检查点。这里只取
    struct_id 最大的那个，也就是跑完整条轨迹之后的终态。
    """
    rows = {}
    with sqlite3.connect(path) as conn:
        for batch_id, batch_name, struct_id, term, value in conn.execute(_SCORE_SQL):
            state = state_of(batch_id, batch_name)
            if state:
                rows.setdefault(state, {}).setdefault(struct_id, {})[term] = value

    out = {}
    for state, by_struct in rows.items():
        out[state] = by_struct[max(by_struct)]        # 最后一个检查点
    return out


def match_state(batch_name):
    """batches.name 形如 ..._bound_wt / ..._unbound_mut，取出其中的状态"""
    name = batch_name.lower()
    for state in ('unbound_wt', 'unbound_mut', 'bound_wt', 'bound_mut'):
        if state in name:              # unbound_* 要排在 bound_* 前面
            return state
    return None


def ddG_from_db3(path):
    """一次轨迹的 ΔΔG -> (GAM 值, 原始值, 逐项 ΔΔG)"""
    states = read_db3(path)
    missing = set(STATE_SIGN) - set(states)
    if missing:
        raise ValueError(f'{path} 缺少状态: {sorted(missing)}；实际有 {sorted(states)}')

    terms = {}
    for state, sign in STATE_SIGN.items():
        for term, value in states[state].items():
            terms[term] = terms.get(term, 0.0) + sign * value

    g, raw = ddG_from_terms(terms)
    return g, raw, terms
