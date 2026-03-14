"""Generate Abaqus vs OpenSees comparison plots.

This standalone script generates overlay plots comparing OpenSees results
with Abaqus reference data for all completed benchmarks.

Usage:
    python scripts/generate_comparison_plots.py

This is a stub — full implementation requires Abaqus reference data files
to be placed in the validation/abaqus_reference/ directory.
"""

from __future__ import annotations

import sys


def main() -> None:
    """Generate comparison plots for all benchmarks."""
    print("Comparison plot generation — not yet implemented.")
    print("Requires Abaqus reference data in validation/abaqus_reference/")
    print("Run Benchmark 01 directly for OpenSees-only result plots:")
    print("  python examples/benchmark_01_abaqus_psi/pipeline_psi_validation.py")


if __name__ == "__main__":
    main()
