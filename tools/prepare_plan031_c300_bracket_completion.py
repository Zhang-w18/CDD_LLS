"""Add auditable post-formal points needed to complete 0.25 dB target brackets."""

from __future__ import annotations

import json
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "outputs" / "experiment031_pdcch_cdd" / "20260907_c300_4tx_2sym"
ADDITIONS = {
    1: {
        "C300_S0_SIDON": [12.25, 12.5, 12.75],
        "C300_AP_RMS_T1": [9.25],
        "C300_AP_TEPS_T1": [12.5],
        "C300_GEO_T1_CTRL_FOLD": [9.25],
        "C300_SMALL_CDD_QSTEP0P25_MATCHED_CDD": [10.5, 10.75, 11.0, 17.25],
        "C300_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD": [11.0],
    },
    2: {
        "C300_B0_QC": [6.0, 6.25],
        "C300_S0_SIDON": [5.75, 6.0],
        "C300_AP_TEPS_T1": [5.25, 5.5],
        "C300_GEO_T1_CTRL_BAL": [6.0, 6.25],
        "C300_GEO_T1_CTRL_FOLD": [2.75],
        "C300_PRG_DFT4_6REG": [9.0],
        "C300_SMALL_CDD_QSTEP0P25_MATCHED_CDD": [10.25],
    },
    4: {
        "C300_AP_RMS_T1": [0.75, 1.0],
        "C300_GEO_T1_CTRL_PAIR": [1.0, 1.25],
        "C300_GEO_T1_CTRL_FOLD": [0.75, 1.0, 1.25, 1.5],
        "C300_PRG_DFT4_6REG": [-1.25],
    },
}


def main() -> None:
    receipt = []
    for al, additions in ADDITIONS.items():
        source = ROOT / "configs" / f"pdcch_result031_c300_al{al}_formal.yaml"
        with open(source, "r", encoding="utf-8") as handle:
            config = yaml.safe_load(handle) or {}
        for candidate in config["candidates"]:
            candidate_id = str(candidate["candidate_id"])
            original = [float(value) for value in candidate["snr_points_db"]]
            added = additions.get(candidate_id, [])
            candidate["snr_points_db"] = sorted(set(original + added))
            if added:
                receipt.append(
                    {
                        "aggregation_level": al,
                        "candidate_id": candidate_id,
                        "added_snr_points_db": added,
                        "reason": "post-formal completion of a <=0.25 dB raw target bracket or bootstrap endpoint support",
                    }
                )
        destination = ROOT / "configs" / f"pdcch_result031_c300_al{al}_bracket_completion.yaml"
        destination.write_text(
            yaml.safe_dump(config, allow_unicode=True, sort_keys=False), encoding="utf-8"
        )
    (OUTPUT / "bracket_completion_receipt.json").write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps({"candidate_shards": len(receipt), "added_points": sum(len(row["added_snr_points_db"]) for row in receipt)}))


if __name__ == "__main__":
    main()
