"""Analyze plan-038 with the requested extended Sidon plotting windows."""

from __future__ import annotations

from run_plan038_refinement import WINDOWS


WINDOWS.update(
    {
        (1, "SIDON_MATCHED"): (4.0, 6.5),
        (2, "SIDON_MATCHED"): (0.0, 2.0),
        (2, "SIDON_TRANSPARENT"): (4.0, 6.0),
    }
)

from analyze_plan038 import main  # noqa: E402


if __name__ == "__main__":
    main()
