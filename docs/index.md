# OpenSees Pipeline Analysis Toolkit

A Python toolkit for nonlinear buried pipeline analysis using OpenSeesPy — grounded in ASCE ALA 2001 and PRCI guidelines, benchmarked against Abaqus.

## What is this?

This toolkit provides a complete workflow for analyzing buried pipelines subjected to geohazard loading. It covers mesh generation along 3D pipeline alignments, ASCE ALA 2001 pipe-soil interaction spring calculations, OpenSees model assembly, analysis execution, and post-processing of results.

The toolkit was developed from professional geotechnical and pipeline engineering practice and is intended as a rigorous, reproducible, open-source alternative to commercial finite element tools like Abaqus for buried pipeline analysis.

## Key Features

- **ALA 2001 soil springs** — Bilinear elastic-perfectly-plastic backbone curves for axial, lateral, bearing, and uplift springs
- **Automated mesh generation** — Node placement along 3D alignments with biased refinement near features
- **Benchmarked against Abaqus** — Validated results within 3.1% of commercial FEA
- **Transparent engineering** — Every equation, assumption, and parameter is documented with references

## Getting Started

```bash
pip install -e ".[dev]"
python examples/benchmark_01_abaqus_psi/pipeline_psi_validation.py
```

See the [Theory](theory/pipe-soil-interaction.md) section for engineering background and the [Validation Report](validation-report.md) for benchmark results.
