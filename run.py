from __future__ import annotations

import argparse
import os

os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/cdd_lls_matplotlib")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

from cdd_lls.core.config import load_config
from cdd_lls.sim.orchestrator import CDDLinkLevelOrchestrator


def main() -> None:
    parser = argparse.ArgumentParser(description="CDD LLS simulation platform")
    parser.add_argument("--config", required=True, help="YAML config path")
    args = parser.parse_args()

    cfg = load_config(args.config)
    runner = CDDLinkLevelOrchestrator(cfg)
    out_dir = runner.run()
    print(f"\nSimulation complete. Results saved to: {out_dir}")


if __name__ == "__main__":
    main()
