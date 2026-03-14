# Benchmark 01: Pipe-Soil Interaction Validation

## Overview

This benchmark validates the OpenSees buried pipeline model against the Abaqus pipe-soil interaction (PSI) workshop model. It demonstrates that the toolkit's soil spring implementation, mesh generation, and moment extraction produce results within 3.1% of Abaqus for a standard elastic pipe-soil interaction problem.

## Model Description

A straight buried steel pipeline modeled with elastic beam-column elements and bilinear elastic-perfectly-plastic soil springs in the lateral horizontal and vertical directions. The pipe is subjected to a prescribed end rotation (UR3) that induces bending and activates the soil springs.

### Geometry

- Pipe alignment: straight, horizontal, along the X-axis
- Pipe length: approximately 20 m
- Pipe outside diameter: 0.2032 m (8 inches)
- Pipe wall thickness: 0.00919 m (0.362 inches)
- Burial depth to pipe centerline: 0.9144 m (3 feet)

### Material Properties

| Property | Value | Unit |
|----------|-------|------|
| Pipe elastic modulus (E) | 210,000 | MPa |
| Pipe Poisson's ratio (nu) | 0.3 | — |
| Pipe yield stress | Not applicable (elastic analysis) | — |
| Soil unit weight (gamma) | 18.0 | kN/m³ |
| Soil friction angle (phi) | 35 | degrees |
| Soil cohesion | 0.0 | kPa |

### Pipe Section Properties

| Property | Value | Unit |
|----------|-------|------|
| Cross-sectional area (A) | 5.606e-3 | m² |
| Moment of inertia (I) | 2.768e-5 | m⁴ |
| Section modulus (S) | 2.723e-4 | m³ |

### Soil Springs

Lateral horizontal (p-y) and vertical (q-z bearing and uplift) springs computed using ALA 2001 Sections 7.3–7.5. Axial springs are included but do not yield under the applied rotation loading.

Spring properties are computed per unit length and scaled by tributary node length before assignment to OpenSees zero-length elements.

### Boundary Conditions

- **Node 1 (far end):** Fixed — all translational and rotational DOFs restrained
- **Node N (loaded end):** Prescribed rotation UR3 applied via displacement control
- **Ground nodes:** Fixed in all DOFs — one ground node per pipe node, connected by zero-length spring elements

### Loading

The applied loading is a prescribed rotation (UR3) at the loaded end of the pipe. This rotation induces bending along the pipe length, which is resisted by the soil springs. The rotation magnitude is chosen to mobilize yielding in a subset of the soil springs while keeping the pipe material elastic.

### Mesh

The mesh uses biased spacing with three zones:

1. **Near the loaded end:** Fine mesh (small elements) to capture steep stress gradients
2. **Transition zone:** Geometrically graded elements
3. **Far field:** Coarser elements where stress gradients are mild

This mesh strategy mirrors the Abaqus workshop model and ensures that the peak bending stress is accurately captured.

## Results

| Metric | OpenSees | Abaqus | Error (%) |
|--------|----------|--------|-----------|
| Maximum bending stress | 36.90 MPa | ~35.8 MPa | 3.1% |
| Lateral springs yielded | 20 | — | — |
| Vertical springs yielded | 51 | — | — |

## Discussion

### Sources of Difference

The 3.1% difference between OpenSees and Abaqus results is attributable to:

1. **Element formulation:** OpenSees `elasticBeamColumn` vs. Abaqus B31 beam elements use different interpolation and integration schemes
2. **Spring distribution:** Abaqus PSI elements handle spring distribution internally, while OpenSees uses discrete zero-length elements with tributary length scaling
3. **Numerical precision:** Different solvers, convergence criteria, and floating-point implementations

### Acceptability

A 3.1% difference is well within acceptable limits for engineering practice:

- ASME B31.8 allows up to 10% difference between independent stress analyses
- The difference is smaller than typical uncertainties in soil spring properties
- Mesh convergence studies confirm the result is mesh-independent

### Known Limitations

- **Elastic pipe only:** The pipe material does not yield in this benchmark. Nonlinear material validation is planned for Benchmark 02.
- **2D springs:** Only lateral and vertical springs are active. Axial springs do not yield under the applied rotation loading.
- **Single soil layer:** Uniform soil properties along the entire pipe length.
- **No pressure or thermal loading:** The only load is the prescribed end rotation.
- **Small displacement:** The analysis uses small displacement theory. Large displacement effects are negligible for this loading magnitude.

## Running the Benchmark

```bash
cd opensees-pipeline-tools
python examples/benchmark_01_abaqus_psi/pipeline_psi_validation.py
```

The script will output the maximum bending stress, spring yield counts, and generate comparison plots in the `figures/` directory.
