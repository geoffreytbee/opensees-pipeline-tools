"""Post-processing results writer for OpenSees pipeline analysis.

This module provides the ``write_results()`` function that standardises output
from any analysis script into a set of CSV and JSON files.  These files can be
loaded by the standalone viewer (``tools/opensees_pipeline_viewer.py``) without
requiring OpenSeesPy.

Output files written by ``write_results()``:

    pipe_nodes.csv
        Nodal coordinates and displacements — one row per pipe node.

    pipe_elements.csv
        Element axial force and end moments — one row per pipe element.

    spring_forces.csv
        Soil spring forces and tributary lengths — one row per node.

    model_metadata.json
        Pipe geometry, material properties, analysis metadata, and any
        additional key/value pairs supplied by the caller.

All spatial coordinates are stored in **metres** regardless of the unit system
used during analysis.  If Imperial units were used, convert to metres *before*
calling ``write_results()``.

This function should be the **last** call in any analysis script (before
``ops.wipe()``).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def write_results(
    output_dir: str,
    x_nodes: np.ndarray,
    disp_ux: np.ndarray,
    disp_uy: np.ndarray,
    disp_uz: np.ndarray,
    elem_My_i: np.ndarray,
    elem_My_j: np.ndarray,
    elem_N: np.ndarray,
    spring_ax: np.ndarray,
    spring_vt: np.ndarray,
    trib: np.ndarray,
    pipe_params: dict[str, float],
    model_metadata: dict[str, Any],
) -> None:
    """Write standardised analysis results for the post-processor viewer.

    This function should be the last call in any analysis script.  It writes
    four output files that the standalone viewer can load without any
    dependency on OpenSeesPy.

    All spatial coordinates must be in metres.  If the analysis used Imperial
    units, convert before calling this function.

    Args:
        output_dir: Directory to write output files.  Created if it does
            not exist.
        x_nodes: Node positions along pipe axis [m], shape ``(n_nodes,)``.
        disp_ux: Nodal displacement in X [m], shape ``(n_nodes,)``.
        disp_uy: Nodal displacement in Y [m], shape ``(n_nodes,)``.
        disp_uz: Nodal displacement in Z [m], shape ``(n_nodes,)``.
        elem_My_i: Bending moment at i-end of each element [N·m],
            shape ``(n_elem,)``.  Sign convention: as returned by
            ``eleForce`` (reaction sign — negate for internal moment).
        elem_My_j: Bending moment at j-end of each element [N·m],
            shape ``(n_elem,)``.
        elem_N: Axial force in each element [N], shape ``(n_elem,)``.
        spring_ax: Axial soil spring force at each node [N],
            shape ``(n_nodes,)``.
        spring_vt: Vertical soil spring force at each node [N],
            shape ``(n_nodes,)``.
        trib: Tributary length at each node [m], shape ``(n_nodes,)``.
        pipe_params: Pipe section properties dictionary with keys:
            ``outer_R`` — outer radius [m],
            ``t_wall`` — wall thickness [m],
            ``E`` — Young's modulus [Pa],
            ``nu`` — Poisson's ratio [-].
        model_metadata: Additional metadata dict.  Common keys include
            ``unit_system``, ``analysis_type``, ``fault_zone_x_min``,
            ``fault_zone_x_max``.  All key/value pairs are written to
            ``model_metadata.json``.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    n_nodes = len(x_nodes)
    n_elem = n_nodes - 1

    # ── pipe_nodes.csv ────────────────────────────────────────────────────
    nodes_df = pd.DataFrame(
        {
            "node_index": np.arange(n_nodes),
            "x_m": x_nodes,
            "UX_m": disp_ux,
            "UY_m": disp_uy,
            "UZ_m": disp_uz,
        }
    )
    nodes_df.to_csv(out / "pipe_nodes.csv", index=False, float_format="%.8e")

    # ── pipe_elements.csv ─────────────────────────────────────────────────
    elems_df = pd.DataFrame(
        {
            "elem_index": np.arange(n_elem),
            "x_i_m": x_nodes[:n_elem],
            "x_j_m": x_nodes[1 : n_elem + 1],
            "N_N": elem_N,
            "My_i_Nm": elem_My_i,
            "My_j_Nm": elem_My_j,
        }
    )
    elems_df.to_csv(out / "pipe_elements.csv", index=False, float_format="%.8e")

    # ── spring_forces.csv ─────────────────────────────────────────────────
    springs_df = pd.DataFrame(
        {
            "node_index": np.arange(n_nodes),
            "x_m": x_nodes,
            "trib_m": trib,
            "axial_spring_N": spring_ax,
            "vertical_spring_N": spring_vt,
        }
    )
    springs_df.to_csv(out / "spring_forces.csv", index=False, float_format="%.8e")

    # ── model_metadata.json ───────────────────────────────────────────────
    outer_r = pipe_params["outer_R"]  # [m]
    t_wall = pipe_params["t_wall"]  # [m]
    e_mod = pipe_params["E"]  # [Pa]
    nu = pipe_params["nu"]  # [-]

    metadata: dict[str, Any] = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "unit_system": model_metadata.get("unit_system", "SI"),
        "pipe_od_m": 2.0 * outer_r,
        "pipe_t_m": t_wall,
        "pipe_E_Pa": e_mod,
        "pipe_nu": nu,
        "analysis_type": model_metadata.get("analysis_type", "unknown"),
    }

    # Merge caller-supplied metadata, converting numpy types for JSON
    for key, value in model_metadata.items():
        if key not in metadata:
            if isinstance(value, np.integer):
                value = int(value)
            elif isinstance(value, np.floating):
                value = float(value)
            elif isinstance(value, np.ndarray):
                value = value.tolist()
            metadata[key] = value

    with open(out / "model_metadata.json", "w") as fh:
        json.dump(metadata, fh, indent=2)
