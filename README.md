# OpenSees Pipeline Analysis Toolkit

**A Python toolkit for nonlinear buried pipeline analysis using OpenSeesPy — grounded in ASCE ALA 2001 and PRCI guidelines, benchmarked against Abaqus**

![Python 3.12+](https://img.shields.io/badge/python-3.12+-blue.svg)
![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)
![CI](https://img.shields.io/github/actions/workflow/status/geoffreybee/opensees-pipeline-tools/ci.yml?label=CI)
![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)

## Overview

OpenSees Pipeline Analysis Toolkit is a Python library for performing nonlinear buried pipeline analysis using [OpenSeesPy](https://openseespydoc.readthedocs.io/). It implements pipe-soil interaction (PSI) models based on ASCE ALA 2001 and PRCI guidelines, providing a rigorous, reproducible open-source alternative to commercial finite element tools for pipeline geohazard analysis. The toolkit grew out of professional geotechnical and pipeline engineering practice and is designed for engineers who need transparent, well-documented computational tools for assessing buried pipeline response to ground movement, fault crossing, slope instability, and seismic loading.

## Engineering Background

Buried pipelines are critical infrastructure that must withstand a range of geohazards including fault crossings, landslides, liquefaction-induced ground movement, and seismic wave propagation. The structural response of a buried pipeline is governed by pipe-soil interaction — the nonlinear force-displacement relationship between the pipe and the surrounding soil. Accurate modeling of PSI requires soil springs in four directions (axial, lateral horizontal, vertical bearing, and uplift) calibrated to site-specific soil conditions. This toolkit implements the industry-standard ASCE ALA 2001 bilinear spring formulations and provides a complete pipeline analysis workflow from mesh generation through post-processing.

## Validation Summary

| Benchmark | Analysis Type | Metric | OpenSees | Abaqus | Error (%) |
|-----------|---------------|--------|----------|--------|-----------|
| Benchmark 01 | Pipe-Soil Interaction | Max bending stress | 36.90 MPa | ~35.8 MPa | 3.1% |

Additional benchmarks covering fault crossing, nonlinear material response, and probabilistic analysis are planned. See the [Validation Report](docs/validation-report.md) for full details.

## Soil Spring Model

The toolkit implements the ASCE ALA 2001 bilinear elastic-perfectly-plastic backbone curves for all four soil spring directions: axial (t-z), lateral horizontal (p-y), vertical bearing (q-z), and uplift (vertical upward). Each spring is characterized by an ultimate force and a yield displacement derived from pipe geometry, burial depth, and soil properties following the equations in ALA 2001 Sections 7.2–7.5. Spring stiffness values are computed per unit length and must be scaled by tributary node length before assignment to OpenSees zero-length elements — this scaling is handled automatically by the toolkit.

## Quick Start

```bash
# Clone the repository
git clone https://github.com/geoffreybee/opensees-pipeline-tools.git
cd opensees-pipeline-tools

# Install in development mode
pip install -e ".[dev]"

# Run Benchmark 01 — pipe-soil interaction validation
python examples/benchmark_01_abaqus_psi/pipeline_psi_validation.py

# Run the test suite
pytest
```

## Project Structure

```
opensees-pipeline-tools/
├── src/opensees_pipeline/      # Core library
│   ├── mesh.py                 # Node/element generation, alignment import
│   ├── springs.py              # ALA 2001 & PRCI soil spring calculations
│   ├── model.py                # Model assembly (nodes, elements, BCs, loads)
│   ├── analysis.py             # Analysis runner, convergence control
│   ├── postprocess.py          # Result extraction (displacements, forces)
│   ├── montecarlo.py           # Monte Carlo framework (planned)
│   └── visualization.py        # 2D/3D result plotting
├── examples/                   # Runnable benchmark models
├── validation/                 # V&V data and comparison results
├── tests/                      # Unit and integration tests
├── docs/                       # MkDocs documentation with theory pages
└── scripts/                    # Standalone utility scripts
```

## Planned Capabilities

- Nonlinear fiber section beam elements (`dispBeamColumn`) for material nonlinearity
- Full ALA 2001 soil spring suite with all four spring directions
- Monte Carlo simulation framework with random field generation
- 3D pipeline visualization with results mapped onto pipe surface (pyvista)
- Shell element models connected to zero-length soil springs
- Seismic wave propagation analysis
- Liquefaction-induced permanent ground displacement (PGD) analysis
- User-selectable unit systems (SI and Imperial)

## References

1. **ASCE ALA (2001).** *Guidelines for the Design of Buried Steel Pipe.* American Lifelines Alliance, American Society of Civil Engineers. — Primary reference for soil spring force-displacement backbone curves (Sections 7.2–7.5).
2. **PRCI.** *Guidelines for the Seismic Design and Assessment of Natural Gas and Liquid Hydrocarbon Pipelines.* Pipeline Research Council International. — Supplementary guidance on pipe-soil interaction and seismic loading.
3. **OpenSeesPy Documentation.** [https://openseespydoc.readthedocs.io/](https://openseespydoc.readthedocs.io/) — Python interface to the OpenSees finite element framework.
4. **Dassault Systèmes (2023).** *Abaqus Analysis User's Manual.* — Reference FEA software used for validation benchmarks.

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

## Author

**Geoffrey Bee** — Geotechnical / Pipeline Engineer

This toolkit was developed from professional engineering practice in buried pipeline analysis and geohazard assessment. Contributions, feedback, and collaboration are welcome.
