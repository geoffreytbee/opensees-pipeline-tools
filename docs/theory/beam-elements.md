# Beam Element Formulations

## Overview

This toolkit supports two beam element formulations for modeling the pipeline:

1. **`elasticBeamColumn`** — linear elastic beam-column elements
2. **`dispBeamColumn`** — displacement-based nonlinear beam-column elements with fiber sections

The choice of element depends on the analysis objective. Linear elastic elements are appropriate for benchmarking, small-displacement problems, and analyses where the pipe remains well within its elastic range. Nonlinear fiber section elements are required when the pipe material yields or when accurate prediction of post-yield behavior is needed.

## elasticBeamColumn

The elastic beam-column element is a standard Euler-Bernoulli (or Timoshenko, if shear deformation is included) beam element with constant section properties. It requires:

- Cross-sectional area $A$
- Moment of inertia $I$ (and $I_z$ for 3D)
- Elastic modulus $E$
- (3D) Shear modulus $G$ and torsional constant $J$

The element produces linear force-displacement relationships and is computationally efficient. It is used in Benchmark 01 for validation against the Abaqus elastic pipe-soil interaction workshop model.

!!! note "Argument Ordering"
    The `elasticBeamColumn` command in OpenSees has specific argument ordering that varies between 2D and 3D models. Always use inline comments to label each argument:
    ```python
    # 2D: tag, iNode, jNode, A, E, I, transfTag
    ops.element('elasticBeamColumn', eleTag, iNode, jNode, A, E, Iz, transfTag)
    ```

## dispBeamColumn with Fiber Sections

The displacement-based beam-column element (`dispBeamColumn`) uses numerical integration along the element length with fiber sections at each integration point. Each fiber has its own uniaxial stress-strain relationship, allowing the element to capture:

- Material nonlinearity (yielding, strain hardening)
- Gradual plastification across the cross-section
- Interaction between axial force and bending moment

For steel pipe sections, the fiber section is typically a circular tube (annular ring) discretized into fibers in the circumferential and radial directions.

This element formulation is planned for future benchmarks involving fault crossing and large ground displacement where pipe yielding is expected.

## When to Use Each Element

| Scenario | Recommended Element |
|----------|-------------------|
| Benchmarking against linear FEA | `elasticBeamColumn` |
| Pipe stress well below yield | `elasticBeamColumn` |
| Fault crossing with large PGD | `dispBeamColumn` + fiber section |
| Slope movement with potential yielding | `dispBeamColumn` + fiber section |
| Thermal + pressure loading near yield | `dispBeamColumn` + fiber section |
| Monte Carlo studies (speed critical) | `elasticBeamColumn` if elastic range confirmed |
