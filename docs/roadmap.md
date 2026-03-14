# Roadmap

## Phase 1: Foundation (Current)

- [x] Repository structure and CI/CD
- [x] ALA 2001 soil spring calculations (all four directions)
- [x] Mesh generation with biased spacing
- [x] Tributary length computation
- [x] Alignment import from CSV/Excel
- [x] Benchmark 01: elastic PSI validation against Abaqus
- [x] Documentation site with MkDocs Material
- [ ] Complete model assembly module (`model.py`)
- [ ] Complete analysis runner (`analysis.py`)
- [ ] Complete post-processing module (`postprocess.py`)

## Phase 2: Nonlinear Material

- [ ] `dispBeamColumn` with fiber section elements
- [ ] Steel material models (bilinear, Ramberg-Osgood)
- [ ] Benchmark 02: nonlinear pipe material validation
- [ ] Combined pressure + bending interaction

## Phase 3: Geohazard Loading

- [ ] Fault crossing analysis with prescribed PGD profile
- [ ] Benchmark 03: fault crossing validation
- [ ] Landslide/slope movement loading
- [ ] Thermal loading (temperature differential)
- [ ] Combined loading scenarios

## Phase 4: Probabilistic Simulation

- [ ] Monte Carlo simulation framework
- [ ] Random field generation (Cholesky decomposition)
- [ ] Random field generation (Karhunen-Loeve expansion)
- [ ] Statistical post-processing and CDF plots
- [ ] Benchmark 04: probabilistic convergence study

## Phase 5: Advanced Modeling

- [ ] PRCI spring model implementation
- [ ] Shell element pipeline models
- [ ] 3D visualization with pyvista
- [ ] Liquefaction-induced PGD analysis
- [ ] Seismic wave propagation analysis
- [ ] User-selectable unit system (SI / Imperial) throughout

## Phase 6: Usability

- [ ] Command-line interface for common workflows
- [ ] Jupyter notebook examples
- [ ] Configuration file (YAML/TOML) for model parameters
- [ ] Automated report generation
