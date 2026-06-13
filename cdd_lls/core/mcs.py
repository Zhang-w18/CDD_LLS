from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List
import math


@dataclass(frozen=True)
class MCS:
    index: int
    qm: int
    code_rate: float
    label: str


# Minimal 38.214-style table entries needed for the QC CDD plan plus a few
# nearby points for debugging. Code rates are normalized by 1024.
NR_256QAM: Dict[int, tuple[int, int]] = {
    0: (2, 120),
    1: (2, 193),
    2: (2, 308),
    3: (2, 449),
    4: (2, 602),
    5: (4, 378),
    6: (4, 434),
    7: (4, 490),
    8: (4, 553),
    9: (4, 616),
    10: (4, 658),
    11: (6, 466),
    12: (6, 517),
    13: (6, 567),
    14: (6, 616),
    15: (6, 666),
    16: (6, 719),
    17: (6, 772),
    18: (6, 822),
    19: (6, 873),
    20: (8, 682),
    21: (8, 711),
    22: (8, 754),
    23: (8, 797),
    24: (8, 841),
    25: (8, 885),
    26: (8, 916),
    27: (8, 948),
}

NR_64QAM: Dict[int, tuple[int, int]] = {
    0: (2, 120),
    1: (2, 193),
    2: (2, 308),
    3: (2, 449),
    4: (2, 602),
    5: (4, 378),
    6: (4, 434),
    7: (4, 490),
    8: (4, 553),
    9: (4, 616),
    10: (4, 658),
    11: (6, 466),
    12: (6, 517),
    13: (6, 567),
    14: (6, 616),
    15: (6, 666),
    16: (6, 719),
    17: (6, 772),
    18: (6, 822),
    19: (6, 873),
    20: (6, 910),
    21: (6, 948),
}


def get_mcs(table: str, index: int, qm: int | None = None, code_rate: float | None = None) -> MCS:
    if qm is not None and code_rate is not None:
        return MCS(index=int(index), qm=int(qm), code_rate=float(code_rate), label="custom")
    name = str(table).lower()
    table_obj = NR_256QAM if name in ("nr_256qam", "256qam") else NR_64QAM
    if int(index) not in table_obj:
        raise ValueError(f"MCS index {index} not available in {table}.")
    qm_val, r1024 = table_obj[int(index)]
    return MCS(index=int(index), qm=int(qm_val), code_rate=float(r1024) / 1024.0, label=name)


@dataclass(frozen=True)
class TBLayout:
    data_re: int
    coded_bits: int
    tb_size: int
    cb_k_values: List[int]
    cb_e_values: List[int]
    qm: int
    code_rate: float


def quantize_tbs(n_info: float) -> int:
    if n_info <= 0:
        return 24
    return max(24, int(8 * round(float(n_info) / 8.0)))


def split_even(total: int, n_parts: int) -> List[int]:
    base = int(total) // int(n_parts)
    rem = int(total) % int(n_parts)
    return [base + (1 if i < rem else 0) for i in range(int(n_parts))]


def build_tb_layout(data_re: int, mcs: MCS, max_cb_payload_bits: int = 3800) -> TBLayout:
    data_re = int(data_re)
    qm = int(mcs.qm)
    coded_bits = int(data_re * qm)
    tb_size = quantize_tbs(coded_bits * float(mcs.code_rate))
    n_cbs = max(1, int(math.ceil(tb_size / float(max_cb_payload_bits))))
    cb_k = split_even(tb_size, n_cbs)
    cb_e = split_even(coded_bits, n_cbs)
    for k, n in zip(cb_k, cb_e):
        if k <= 0 or n <= 0 or k >= n:
            raise ValueError(f"Invalid CB sizing: k={k}, n={n}, coded_bits={coded_bits}, tb_size={tb_size}")
        rate = float(k) / float(n)
        if rate < 0.05 or rate > 0.95:
            raise ValueError(f"Unsupported toy LDPC rate {rate:.3f} for k={k}, n={n}")
    return TBLayout(
        data_re=data_re,
        coded_bits=coded_bits,
        tb_size=tb_size,
        cb_k_values=cb_k,
        cb_e_values=cb_e,
        qm=qm,
        code_rate=float(mcs.code_rate),
    )
