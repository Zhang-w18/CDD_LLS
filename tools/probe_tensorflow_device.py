"""Print a JSON receipt for TensorFlow CPU/GPU visibility and basic placement."""

from __future__ import annotations

import json

from tools.export_bler_result_bundle import environment_receipt


if __name__ == "__main__":
    print(json.dumps(environment_receipt(), indent=2, ensure_ascii=False))
