# Expected Output — Benchmark 01

## Key Results

| Metric | Expected Value |
|--------|---------------|
| Maximum bending stress | 36.90 MPa |
| Location of max stress | Near the loaded end (~0.3 m) |
| Lateral springs yielded | ~20 out of 34 |
| Vertical springs yielded | ~51 out of 68 (bearing + uplift) |
| Analysis convergence | Successful in 100 load steps |

## Validation Criteria

The benchmark is considered passing if:

1. Maximum bending stress is within 5% of 36.90 MPa
2. The analysis converges without warnings
3. The moment diagram shape is physically reasonable (peak near loaded end, monotonically decaying)
4. Sum of tributary lengths equals total pipe length (within floating-point tolerance)

## Comparison with Abaqus

| Metric | OpenSees | Abaqus | Error |
|--------|----------|--------|-------|
| Max bending stress | 36.90 MPa | ~35.8 MPa | 3.1% |

The 3.1% difference is within acceptable engineering tolerance and is attributed to differences in element formulation, integration schemes, and spring distribution methods between the two codes.
