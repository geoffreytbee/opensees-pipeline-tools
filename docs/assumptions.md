# Engineering Assumptions Register

This page documents the key engineering assumptions and simplifications made in the toolkit. Understanding these assumptions is essential for interpreting results and determining the scope of validity.

## Soil Spring Assumptions

| # | Assumption | Justification | Impact if Violated |
|---|-----------|---------------|-------------------|
| A1 | Soil springs are independent in each direction | ALA 2001 convention; coupling effects are second-order for typical buried pipeline problems | May underestimate resistance for combined loading at large displacements |
| A2 | Bilinear elastic-perfectly-plastic backbone | ALA 2001 standard model; captures essential nonlinearity with minimum parameters | Does not capture strain softening, cyclic degradation, or rate effects |
| A3 | Spring properties are uniform along the pipe unless varied by user | Simplification for initial implementation | Must be manually varied for layered soil profiles |
| A4 | No gap formation behind the pipe | Assumed soil remains in contact with pipe on all sides | Conservative for uplift; unconservative for lateral if gap forms |
| A5 | Springs act at pipe centerline | Standard practice; acceptable when pipe diameter is small relative to burial depth | May introduce error for very shallow or large-diameter pipes |

## Structural Assumptions

| # | Assumption | Justification | Impact if Violated |
|---|-----------|---------------|-------------------|
| S1 | Pipe behaves as a beam (plane sections remain plane) | Valid for D/t > 20 and no local buckling | Shell elements required for thick-walled pipes or local buckling assessment |
| S2 | Small strain theory for elastic analyses | Benchmark 01 displacements are small relative to pipe length | Large displacement formulation needed for fault crossing with large PGD |
| S3 | No initial imperfections or residual stresses | Simplification for initial validation | May affect buckling-critical analyses |
| S4 | Pipe cross-section remains circular | No ovalization modeled with beam elements | Shell elements required for ovalization assessment |

## Loading Assumptions

| # | Assumption | Justification | Impact if Violated |
|---|-----------|---------------|-------------------|
| L1 | Quasi-static loading | Inertial effects neglected; valid for slow ground movement | Dynamic analysis required for seismic wave propagation |
| L2 | Monotonic loading only | No load reversal or cyclic degradation | Cyclic material models needed for seismic analysis |
| L3 | No internal pressure effects on spring stiffness | Pressure does not significantly affect soil spring properties | Pressure affects pipe stiffness and may change failure mode |

## Numerical Assumptions

| # | Assumption | Justification | Impact if Violated |
|---|-----------|---------------|-------------------|
| N1 | Convergence tolerance of 1e-6 on norm of unbalanced force | Standard practice for nonlinear FEA | Tighter tolerance may be needed for highly nonlinear problems |
| N2 | Tributary length scaling produces equivalent continuous spring model | Mathematically exact as mesh is refined | Coarse meshes may introduce discretization error in spring forces |
