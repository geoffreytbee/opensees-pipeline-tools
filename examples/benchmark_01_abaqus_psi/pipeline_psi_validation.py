"""Benchmark 01: Pipe-Soil Interaction Validation against Abaqus.

This script builds and runs an OpenSees model of a straight buried pipeline
with elastic beam-column elements and ALA 2001 bilinear soil springs,
replicating the Abaqus PSI workshop model.

Validated result: maximum bending stress = 36.90 MPa (Abaqus: ~35.8 MPa, error: 3.1%)

Usage:
    python examples/benchmark_01_abaqus_psi/pipeline_psi_validation.py
"""

from __future__ import annotations

import math
import sys

import numpy as np

try:
    import openseespy.opensees as ops
except ImportError:
    print("ERROR: openseespy is required. Install with: pip install openseespy")
    sys.exit(1)

import matplotlib.pyplot as plt

# Add src to path for library imports
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[2] / "src"))

from opensees_pipeline.mesh import biased_spacing, compute_tributary_lengths
from opensees_pipeline.springs import ala2001_lateral_spring, ala2001_bearing_spring, ala2001_uplift_spring


def run_benchmark_01() -> None:
    """Run the Benchmark 01 pipe-soil interaction validation model."""

    # =========================================================================
    # MODEL PARAMETERS
    # =========================================================================

    # Pipe geometry
    pipe_od = 0.2032  # outside diameter [m] (8 inches)
    pipe_wt = 0.00919  # wall thickness [m] (0.362 inches)
    pipe_id = pipe_od - 2.0 * pipe_wt  # inside diameter [m]

    # Pipe section properties
    area = math.pi / 4.0 * (pipe_od**2 - pipe_id**2)  # cross-sectional area [m^2]
    inertia = math.pi / 64.0 * (pipe_od**4 - pipe_id**4)  # moment of inertia [m^4]

    # Pipe material (elastic)
    e_steel = 210_000_000.0  # elastic modulus [kPa] (210 GPa)

    # Soil properties
    burial_depth = 0.9144  # depth to pipe centerline [m] (3 ft)
    soil_gamma = 18.0  # soil unit weight [kN/m^3]
    phi = 35.0  # friction angle [degrees]

    # =========================================================================
    # MESH GENERATION
    # =========================================================================

    # Three mesh zones mirroring Abaqus workshop model
    # Zone 1: fine mesh near loaded end (0 to 5 m)
    zone1 = biased_spacing(n_elem=15, length=5.0, bias_ratio=2.0, fine_end="start")
    # Zone 2: transition (5 to 10 m)
    zone2 = biased_spacing(n_elem=10, length=5.0, bias_ratio=1.5, fine_end="start") + 5.0
    # Zone 3: coarse far field (10 to 20 m)
    zone3 = biased_spacing(n_elem=10, length=10.0, bias_ratio=1.0) + 10.0

    # Combine zones, removing duplicate nodes at zone boundaries
    x_nodes = np.concatenate([zone1, zone2[1:], zone3[1:]])
    n_nodes = len(x_nodes)
    n_elem = n_nodes - 1

    print(f"Mesh: {n_nodes} nodes, {n_elem} elements")
    print(f"Pipe length: {x_nodes[-1]:.2f} m")

    # Tributary lengths for spring scaling
    trib_lengths = compute_tributary_lengths(x_nodes)
    print(f"Sum of tributary lengths: {trib_lengths.sum():.4f} m (should equal pipe length)")

    # =========================================================================
    # SOIL SPRING PROPERTIES (per unit length)
    # =========================================================================

    pu, yu = ala2001_lateral_spring(pipe_od, burial_depth, soil_gamma, phi)
    qu_bear, zq_bear = ala2001_bearing_spring(pipe_od, burial_depth, soil_gamma, phi)
    qu_uplift, zq_uplift = ala2001_uplift_spring(pipe_od, burial_depth, soil_gamma, phi)

    print(f"\nSoil spring properties (per unit length):")
    print(f"  Lateral:  pu = {pu:.2f} kN/m,  yu = {yu:.4f} m")
    print(f"  Bearing:  qu = {qu_bear:.2f} kN/m,  zq = {zq_bear:.4f} m")
    print(f"  Uplift:   qu = {qu_uplift:.2f} kN/m,  zq = {zq_uplift:.4f} m")

    # =========================================================================
    # OPENSEES MODEL
    # =========================================================================

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)  # 2D model, 3 DOFs per node (ux, uy, rz)

    # --- Nodes ---
    # Pipe nodes: tags 1 to n_nodes
    for i in range(n_nodes):
        ops.node(
            i + 1,  # nodeTag
            x_nodes[i],  # x coordinate [m]
            0.0,  # y coordinate [m]
        )

    # Ground nodes for springs: tags n_nodes+1 to 2*n_nodes
    for i in range(n_nodes):
        ops.node(
            n_nodes + i + 1,  # nodeTag
            x_nodes[i],  # x coordinate [m]
            0.0,  # y coordinate [m]
        )
        ops.fix(n_nodes + i + 1, 1, 1, 1)  # fix all DOFs on ground nodes

    # --- Boundary Conditions ---
    # Fixed end (far end, last node)
    ops.fix(n_nodes, 1, 1, 1)  # fix ux, uy, rz

    # --- Materials for soil springs ---
    # Material tags: lateral = 100+i, bearing = 200+i, uplift = 300+i
    for i in range(n_nodes):
        trib = trib_lengths[i]  # tributary length [m]

        # Lateral spring (horizontal, DOF 1 = ux)
        pu_node = pu * trib  # ultimate force at this node [kN]
        ky_lat = pu_node / yu  # initial stiffness [kN/m]
        ops.uniaxialMaterial(
            "ElasticPP",  # bilinear elastic-perfectly-plastic
            100 + i + 1,  # matTag
            ky_lat,  # E: initial stiffness [kN/m]
            yu,  # epsyP: yield displacement [m]
        )

        # Bearing spring (vertical downward, positive y direction mapped to DOF 2)
        qu_bear_node = qu_bear * trib  # ultimate bearing force [kN]
        kz_bear = qu_bear_node / zq_bear  # initial stiffness [kN/m]
        ops.uniaxialMaterial(
            "ElasticPP",
            200 + i + 1,  # matTag
            kz_bear,  # E [kN/m]
            zq_bear,  # epsyP [m]
        )

        # Uplift spring (vertical upward)
        qu_up_node = qu_uplift * trib  # ultimate uplift force [kN]
        kz_up = qu_up_node / zq_uplift  # initial stiffness [kN/m]
        ops.uniaxialMaterial(
            "ElasticPP",
            300 + i + 1,  # matTag
            kz_up,  # E [kN/m]
            zq_uplift,  # epsyP [m]
        )

    # --- Zero-length spring elements ---
    # Element tags: springs start at 1000
    spring_tag = 1000
    for i in range(n_nodes):
        pipe_node = i + 1
        ground_node = n_nodes + i + 1

        # Lateral spring (DOF 1 = ux)
        ops.element(
            "zeroLength",
            spring_tag,  # eleTag
            ground_node,  # iNode (ground)
            pipe_node,  # jNode (pipe)
            "-mat", 100 + i + 1,  # material tag
            "-dir", 1,  # DOF direction (1 = ux)
        )
        spring_tag += 1

        # Vertical spring (DOF 2 = uy) — combined bearing/uplift via bearing material
        ops.element(
            "zeroLength",
            spring_tag,  # eleTag
            ground_node,  # iNode (ground)
            pipe_node,  # jNode (pipe)
            "-mat", 200 + i + 1,  # material tag
            "-dir", 2,  # DOF direction (2 = uy)
        )
        spring_tag += 1

    # --- Geometric transformation ---
    ops.geomTransf("Linear", 1)  # linear geometric transformation, tag=1

    # --- Pipe beam-column elements ---
    for i in range(n_elem):
        ops.element(
            "elasticBeamColumn",
            i + 1,  # eleTag
            i + 1,  # iNode
            i + 2,  # jNode
            area,  # A [m^2]
            e_steel,  # E [kPa]
            inertia,  # I [m^4]
            1,  # transfTag
        )

    # =========================================================================
    # ANALYSIS — DISPLACEMENT CONTROL
    # =========================================================================

    # Apply prescribed rotation at node 1 (loaded end)
    target_rotation = 0.01  # target rotation [rad] at loaded end
    n_steps = 100

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.sp(1, 3, target_rotation)  # node 1, DOF 3 (rz), target value

    ops.constraints("Penalty", 1.0e14, 1.0e14)
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormUnbalance", 1.0e-6, 100)  # tolerance, max iterations
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / n_steps)  # delta lambda
    ops.analysis("Static")

    print(f"\nRunning analysis: {n_steps} load steps...")
    ok = ops.analyze(n_steps)

    if ok != 0:
        print("WARNING: Analysis did not converge!")
    else:
        print("Analysis converged successfully.")

    # =========================================================================
    # POST-PROCESSING
    # =========================================================================

    # Extract bending moments using eleForce (negate My_i)
    moments = np.zeros(n_elem)
    for i in range(n_elem):
        # eleForce returns [N_i, V_i, M_i, N_j, V_j, M_j] for 2D beam
        forces = ops.eleForce(i + 1)
        moments[i] = -forces[2]  # negate M_i [kN*m]

    # Bending stress: sigma = M * (D/2) / I
    stresses = np.abs(moments) * (pipe_od / 2.0) / inertia  # [kPa]
    stresses_mpa = stresses / 1000.0  # convert kPa to MPa

    max_stress = np.max(stresses_mpa)
    max_stress_loc = x_nodes[np.argmax(stresses_mpa)]

    print(f"\n{'='*60}")
    print(f"RESULTS")
    print(f"{'='*60}")
    print(f"Maximum bending stress: {max_stress:.2f} MPa")
    print(f"Location of max stress: x = {max_stress_loc:.2f} m")
    print(f"Abaqus reference:       ~35.8 MPa")
    print(f"Error:                  {abs(max_stress - 35.8) / 35.8 * 100:.1f}%")

    # Extract displacements
    disp_x = np.zeros(n_nodes)
    disp_y = np.zeros(n_nodes)
    for i in range(n_nodes):
        disp_x[i] = ops.nodeDisp(i + 1, 1)  # ux [m]
        disp_y[i] = ops.nodeDisp(i + 1, 2)  # uy [m]

    # Count yielded springs
    n_lat_yielded = 0
    n_vert_yielded = 0
    for i in range(n_nodes):
        # Check if lateral spring has yielded (displacement > yu)
        lat_disp = abs(ops.nodeDisp(i + 1, 1))
        if lat_disp > yu:
            n_lat_yielded += 1
        # Check vertical
        vert_disp = abs(ops.nodeDisp(i + 1, 2))
        if vert_disp > min(zq_bear, zq_uplift):
            n_vert_yielded += 1

    print(f"Lateral springs yielded:  {n_lat_yielded} / {n_nodes}")
    print(f"Vertical springs yielded: {n_vert_yielded} / {n_nodes}")

    # =========================================================================
    # PLOTTING
    # =========================================================================

    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)

    # Displacement profile
    axes[0].plot(x_nodes, disp_y * 1000, "b-", linewidth=1.5)
    axes[0].set_ylabel("Vertical displacement [mm]")
    axes[0].set_title("Benchmark 01: Pipe-Soil Interaction Validation")
    axes[0].grid(True, alpha=0.3)

    # Bending moment diagram
    axes[1].plot(x_nodes[:n_elem], moments, "r-", linewidth=1.5)
    axes[1].set_ylabel("Bending moment [kN*m]")
    axes[1].grid(True, alpha=0.3)

    # Bending stress profile
    axes[2].plot(x_nodes[:n_elem], stresses_mpa, "g-", linewidth=1.5)
    axes[2].axhline(y=35.8, color="k", linestyle="--", alpha=0.5, label="Abaqus (~35.8 MPa)")
    axes[2].set_ylabel("Bending stress [MPa]")
    axes[2].set_xlabel("Distance along pipe [m]")
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()

    # Save figure
    import pathlib
    fig_dir = pathlib.Path(__file__).resolve().parents[2] / "figures"
    fig_dir.mkdir(exist_ok=True)
    fig_path = fig_dir / "benchmark_01_results.png"
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    print(f"\nPlot saved to: {fig_path}")
    plt.close()

    ops.wipe()


if __name__ == "__main__":
    run_benchmark_01()
