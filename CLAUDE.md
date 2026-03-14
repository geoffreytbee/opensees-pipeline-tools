# CLAUDE.md — Project Memory for Claude Code

## Project Summary

OpenSees Pipeline Analysis Toolkit: a Python library for nonlinear buried pipeline
analysis using OpenSeesPy. Focus areas are pipe-soil interaction (PSI), geohazard
loading (fault crossing, landslide, seismic, liquefaction PGD), and probabilistic
simulation. Benchmarked against Abaqus. Grounded in ASCE ALA 2001 and PRCI guidelines.

## Stack

- Python 3.12
- openseespy — finite element analysis engine
- numpy — numerical computation
- matplotlib — 2D plotting
- pandas — data handling and CSV/Excel I/O
- openpyxl — Excel file support
- pyvista — 3D visualization (planned, not yet integrated)
- mkdocs-material + mkdocstrings — documentation site

## Environment

- Virtual environment: `~/Documents/Claude-Projects/os_env`
- Activate before working: `source ~/Documents/Claude-Projects/os_env/bin/activate`
- Install in dev mode: `pip install -e ".[dev]"`

## Git Conventions

- Conventional commits: `feat:`, `fix:`, `docs:`, `test:`, `refactor:`, `chore:`
- Feature branches off `main`; never commit directly to `main`
- PRs require passing CI (ruff, black, pytest)

## Code Style

- Google-style docstrings on all public functions and classes
- Type hints on all function signatures
- `snake_case` for functions, variables, modules; `PascalCase` for classes
- No wildcard imports (`from x import *`)
- No bare `except:` clauses
- Line length: 88 (black default)
- Linter: ruff (E, F, I, W rules)

## Critical Engineering Conventions

### Unit Labeling
- Every numerical constant must have its unit in a comment
- SI default (m, kN, kPa, MPa) unless user selects Imperial
- All functions that accept dimensional quantities must document units in docstring

### Spring Stiffness Convention
- **All spring stiffness/force values passed to OpenSees are PER NODE** (scaled by
  tributary length), NEVER per unit length
- The `springs.py` functions return per-unit-length values
- Scaling by tributary length happens during model assembly in `model.py`
- This is a critical distinction — getting it wrong silently produces wrong results

### Element Convention
- `elasticBeamColumn` — for linear/benchmark analyses (Benchmark 01 uses this)
- `dispBeamColumn` + fiber section — for nonlinear material analyses
- Both element commands have non-obvious argument ordering; always use inline comments
  labeling each argument

### Moment Extraction Convention
- Use `eleForce` with node-based extraction
- Negate `My_i` to get correct sign for bending moment at node i
- **DO NOT** average `My_j[left]` and `My_i[right]` at shared nodes — this causes
  cancellation and produces near-zero values at interior nodes
- Extract from a single element per node (use element to the right of each node)

## Validated Benchmarks

| Benchmark | Metric | OpenSees | Abaqus | Error |
|-----------|--------|----------|--------|-------|
| 01 — PSI | Max bending stress | 36.90 MPa | ~35.8 MPa | 3.1% |

## File Organization

- `src/opensees_pipeline/` — library code (mesh, springs, model, analysis, postprocess)
- `examples/` — runnable benchmark scripts with READMEs
- `validation/` — V&V data, Abaqus reference files, comparison results
- `tests/` — unit tests mirroring `src/` structure
- `docs/` — MkDocs pages (theory, architecture, validation report)
- `scripts/` — standalone utility scripts
- `figures/` — generated plots and figures

## Common Commands

```bash
# Run tests
pytest

# Lint
ruff check .

# Format
black .

# Check formatting without modifying
black --check .

# Serve docs locally
mkdocs serve

# Deploy docs to GitHub Pages
mkdocs gh-deploy --force

# Run Benchmark 01
python examples/benchmark_01_abaqus_psi/pipeline_psi_validation.py
```
