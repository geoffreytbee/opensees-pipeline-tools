# Verification and Validation (V&V) Methodology

## Purpose

This document describes the verification and validation methodology used to establish confidence in the OpenSees Pipeline Analysis Toolkit. V&V is essential for computational engineering tools because errors in the code or modeling approach can lead to unsafe design decisions.

## Verification vs. Validation

**Verification** answers: "Are we solving the equations correctly?"

- Checks that the code correctly implements the intended mathematical model
- Compares against analytical solutions, manufactured solutions, or mesh convergence studies
- Focuses on numerical accuracy, convergence, and implementation correctness

**Validation** answers: "Are we solving the right equations?"

- Checks that the mathematical model adequately represents the physical problem
- Compares against experimental data, field measurements, or trusted reference codes
- Focuses on the appropriateness of assumptions, idealizations, and constitutive models

## Benchmark Structure

Each benchmark follows a standardized structure:

### 1. Problem Description

- Geometry: pipe alignment, diameter, wall thickness, burial depth
- Materials: pipe material properties, soil properties
- Loading: type (displacement, force, thermal), magnitude, application method
- Boundary conditions: fixity, symmetry, spring supports

### 2. Engineering Assumptions

- Beam vs. shell idealization
- Small vs. large displacement theory
- Constitutive model selection (elastic, bilinear, fiber section)
- Soil spring model (ALA 2001, PRCI, or custom)
- Drainage condition (drained vs. undrained soil)

### 3. Mesh and Element Formulation

- Element type and formulation
- Mesh density and refinement strategy
- Mesh convergence evidence
- Integration scheme (number of integration points)

### 4. Solver Configuration

- Analysis type (static, dynamic, eigenvalue)
- Constraint handler
- Numberer
- System of equations solver
- Convergence test (norm type, tolerance, maximum iterations)
- Solution algorithm (Newton, Modified Newton, Broyden)
- Integrator (load control, displacement control, arc-length)

### 5. Reference Solution

- Source: analytical, Abaqus, experimental, published literature
- Key reference values with units
- Conditions under which the reference solution is valid

### 6. Comparison Methodology

- Quantities compared (stress, displacement, reaction force, etc.)
- Error metric: relative error = |OpenSees - Reference| / |Reference| x 100%
- Acceptance criteria: typically < 5% for engineering validation, < 1% for verification against analytical solutions

### 7. Results and Discussion

- Quantitative comparison table
- Plots overlaying OpenSees and reference results
- Discussion of sources of difference
- Statement on scope of validity

### 8. Known Limitations

- What the benchmark does NOT validate
- Conditions under which the validated model should not be used
- Recommended follow-up benchmarks

## Current Benchmarks

| ID | Description | Status | Max Error |
|----|-------------|--------|-----------|
| 01 | Elastic PSI — Abaqus workshop model | Complete | 3.1% |
| 02 | Nonlinear pipe — fault crossing | Planned | — |
| 03 | Probabilistic convergence study | Planned | — |

## Directory Structure

```
validation/
├── README.md               # This file
├── abaqus_reference/       # Abaqus reference data (.odb excluded from git)
│   └── .gitkeep
└── comparison/             # Comparison scripts and results
    └── .gitkeep
```

Reference data from Abaqus (ODB files, extracted CSV data) should be placed in `abaqus_reference/`. Comparison scripts that generate overlay plots and error tables go in `comparison/`.
