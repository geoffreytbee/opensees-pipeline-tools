# Abaqus Reference Model

This directory is reserved for reference data from the Abaqus PSI workshop model used as the validation target for Benchmark 01.

## Abaqus Workshop Model

The reference model is based on the Abaqus pipe-soil interaction workshop example, which demonstrates the use of PSI (Pipe-Soil Interaction) elements for buried pipeline analysis. Key features:

- **Element type:** B31 (Timoshenko beam) for pipe, PSI (pipe-soil interaction) for soil springs
- **Material:** Linear elastic steel
- **Loading:** Prescribed rotation (UR3) at one end
- **Boundary conditions:** Fixed at the opposite end
- **Soil springs:** Bilinear elastic-perfectly-plastic, lateral and vertical directions

## Reference Result

Maximum bending stress: approximately 35.8 MPa

## File Placement

Abaqus output database (.odb) files and result extractions should be placed in this directory. These files are excluded from version control via `.gitignore` due to their binary format and large size.
