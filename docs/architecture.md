# Architecture

## Project Structure

```
src/opensees_pipeline/
├── mesh.py           # Mesh generation and alignment import
├── springs.py        # ALA 2001 soil spring backbone curves
├── model.py          # Model assembly (nodes, elements, BCs, loads)
├── analysis.py       # Analysis execution and convergence control
├── postprocess.py    # Result extraction and stress computation
├── montecarlo.py     # Monte Carlo framework (planned)
└── visualization.py  # 2D/3D result visualization
```

## Data Flow

A typical analysis follows this pipeline:

```
Alignment CSV/Excel
        │
        ▼
    mesh.py ──── read_alignment_csv(), generate_mesh_from_kps()
        │         compute_tributary_lengths()
        ▼
  springs.py ──── ala2001_axial_spring(), ala2001_lateral_spring(), etc.
        │         Returns per-unit-length values
        ▼
   model.py ──── Builds OpenSees model: nodes, elements, BCs
        │         Scales springs by tributary length (per-unit-length → per-node)
        │         Attaches zero-length spring elements
        ▼
 analysis.py ──── Configures solver, runs load steps
        │         Handles convergence with adaptive sub-stepping
        ▼
postprocess.py ── Extracts displacements, element forces, stresses
        │
        ▼
visualization.py ── Plots displacement profiles, moment diagrams, etc.
```

## Key Design Decisions

1. **Per-unit-length springs in `springs.py`**: Spring functions return values per unit length, not per node. This keeps the spring calculations independent of the mesh and makes them reusable. Tributary length scaling is performed in `model.py` during assembly.

2. **Separation of mesh and model**: Mesh generation (`mesh.py`) is independent of OpenSees. It produces numpy arrays of coordinates and connectivity that can be used by any FE code. The OpenSees-specific model construction lives in `model.py`.

3. **Moment extraction convention**: Bending moments are extracted per-element at node i (negated) to avoid cancellation at shared interior nodes. This is enforced in `postprocess.py`.

4. **Unit system as parameter**: Functions accept a `unit_system` parameter rather than performing internal conversions. The user is responsible for providing consistent units.
