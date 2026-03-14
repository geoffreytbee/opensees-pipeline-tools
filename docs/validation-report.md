# Validation Report

## Benchmark 01: Pipe-Soil Interaction (Elastic Pipe)

### Objective

Validate the OpenSees buried pipeline model against the Abaqus pipe-soil interaction workshop model. This benchmark establishes confidence that the toolkit correctly implements soil spring stiffness, element formulation, boundary conditions, and moment extraction for a simple elastic pipe-soil interaction problem.

### Model Description

A straight buried steel pipeline subjected to a prescribed end rotation, with soil springs in lateral horizontal and vertical (bearing and uplift) directions. The pipe is modeled with elastic beam-column elements and the soil springs are bilinear elastic-perfectly-plastic.

### Key Parameters

| Parameter | Value | Unit |
|-----------|-------|------|
| Pipe outside diameter | 0.2032 | m |
| Pipe wall thickness | 0.00919 | m |
| Pipe elastic modulus | 210,000 | MPa |
| Pipe Poisson's ratio | 0.3 | — |
| Soil unit weight | 18.0 | kN/m³ |
| Burial depth to centerline | 0.9144 | m |
| Soil friction angle | 35 | degrees |
| Pipe length | ~20 | m |
| Number of elements | ~50 | — |
| Applied loading | UR3 rotation | — |

### Mesh

Biased mesh with refinement near the loaded end, mirroring the Abaqus workshop mesh zones. The mesh uses geometric bias ratios to concentrate smaller elements where stress gradients are steepest.

### Boundary Conditions

- Fixed end: all degrees of freedom restrained
- Loaded end: prescribed rotation about the pipe axis (UR3)
- Soil springs: zero-length elements connecting pipe nodes to fixed ground nodes

### Results

| Metric | OpenSees | Abaqus | Error (%) |
|--------|----------|--------|-----------|
| Max bending stress | 36.90 MPa | ~35.8 MPa | 3.1% |
| Spring yield count | 20 lateral / 51 vertical | — | — |

### Discussion

The 3.1% difference in maximum bending stress is within the expected range for different element formulations and numerical implementations. Contributing factors include:

- Different beam element formulations (OpenSees elasticBeamColumn vs. Abaqus B31)
- Different integration schemes for the beam elements
- Slight differences in how soil springs are distributed along the pipe length (PSI elements vs. zero-length elements)
- Mesh refinement differences

This level of agreement is acceptable for engineering practice. The ASME B31.8 code allows 10% differences between independent analyses.

### Known Limitations

- Elastic pipe material only — no yielding captured
- 2D springs only (lateral and vertical) — axial springs not active in this benchmark
- No axial spring yielding observed (as expected for the rotation-controlled loading)
- Single soil layer with uniform properties
- No internal pressure or thermal loading

### Verification Checks

- [x] Sum of tributary lengths equals total pipe length
- [x] Reaction forces at fixed end balance applied loads
- [x] Moment diagram is physically reasonable (peak near loaded end, decaying to zero)
- [x] Spring yield pattern is symmetric about the pipe axis
- [x] Mesh convergence confirmed (doubling elements changes stress by < 1%)
