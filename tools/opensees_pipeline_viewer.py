#!/usr/bin/env python3
"""OpenSees Pipeline Viewer — standalone 3D post-processor.

Loads results written by ``opensees_pipeline.postprocess.write_results()`` and
renders an interactive 3D visualisation of the pipeline with result contours
plus a 2D companion profile plot.

Dependencies: numpy, matplotlib, pyvista, pandas
Does **not** require openseespy.

Usage::

    python tools/opensees_pipeline_viewer.py [results_dir]

If *results_dir* is not provided, the script prompts for drag-and-drop input.
"""

from __future__ import annotations

import json
import os
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import pyvista as pv
except ImportError:
    print("ERROR: pyvista is not installed.")
    print("Install it with: pip install pyvista")
    sys.exit(1)

import matplotlib.pyplot as plt

# ── Constants ─────────────────────────────────────────────────────────────

N_CIRC = 16  # circumferential discretisation points per ring

RESULT_INFO: dict[str, dict[str, str]] = {
    "sigma_long": {"label": "Longitudinal Stress", "unit": "Pa", "key": "1"},
    "sigma_vm": {"label": "Von Mises Stress", "unit": "Pa", "key": "2"},
    "sigma_hoop": {"label": "Hoop Stress", "unit": "Pa", "key": "3"},
    "epsilon_long": {"label": "Longitudinal Strain", "unit": "m/m", "key": "4"},
    "disp_vert": {"label": "Vertical Displacement", "unit": "m", "key": "5"},
    "disp_ax": {"label": "Axial Displacement", "unit": "m", "key": "6"},
    "moment": {"label": "Bending Moment", "unit": "N·m", "key": "7"},
    "spring_vt": {"label": "Vertical Spring Force", "unit": "N", "key": "8"},
    "spring_ax": {"label": "Axial Spring Force", "unit": "N", "key": "9"},
}

RESULT_KEYS = list(RESULT_INFO.keys())

# Display-unit conversion factors (raw SI → friendly display units)
_DISPLAY_SCALE: dict[str, float] = {
    "sigma_long": 1e-6,  # Pa → MPa
    "sigma_vm": 1e-6,  # Pa → MPa
    "sigma_hoop": 1e-6,  # Pa → MPa
    "epsilon_long": 1e3,  # m/m → mm/m (millistrain)
    "disp_vert": 1e3,  # m → mm
    "disp_ax": 1e3,  # m → mm
    "moment": 1e-3,  # N·m → kN·m
    "spring_vt": 1e-3,  # N → kN
    "spring_ax": 1e-3,  # N → kN
}

_DISPLAY_YLABEL: dict[str, str] = {
    "sigma_long": "Stress [MPa]",
    "sigma_vm": "Stress [MPa]",
    "sigma_hoop": "Stress [MPa]",
    "epsilon_long": "Strain [mm/m]",
    "disp_vert": "Displacement [mm]",
    "disp_ax": "Displacement [mm]",
    "moment": "Moment [kN·m]",
    "spring_vt": "Force [kN]",
    "spring_ax": "Force [kN]",
}


# ── Data Loading ──────────────────────────────────────────────────────────


def load_results(results_dir: str) -> dict:
    """Load all output files from a results directory.

    Args:
        results_dir: Path to directory containing ``pipe_nodes.csv``,
            ``pipe_elements.csv``, ``spring_forces.csv``, and
            ``model_metadata.json``.

    Returns:
        Dictionary with keys ``nodes_df``, ``elems_df``, ``springs_df``,
        ``metadata``, and ``results_dir``.

    Raises:
        FileNotFoundError: If any required file is missing.
    """
    rd = Path(results_dir)

    required = [
        "pipe_nodes.csv",
        "pipe_elements.csv",
        "spring_forces.csv",
        "model_metadata.json",
    ]
    for fname in required:
        if not (rd / fname).exists():
            raise FileNotFoundError(f"Missing required file: {rd / fname}")

    nodes_df = pd.read_csv(rd / "pipe_nodes.csv")
    elems_df = pd.read_csv(rd / "pipe_elements.csv")
    springs_df = pd.read_csv(rd / "spring_forces.csv")

    with open(rd / "model_metadata.json") as fh:
        metadata = json.load(fh)

    return {
        "nodes_df": nodes_df,
        "elems_df": elems_df,
        "springs_df": springs_df,
        "metadata": metadata,
        "results_dir": str(rd),
    }


# ── Derived Quantities ────────────────────────────────────────────────────


def compute_derived(data: dict) -> dict:
    """Compute derived stress, strain, and section quantities.

    Populates ``data`` with additional keys used by the surface builder
    and scalar field calculations.

    Args:
        data: Dictionary returned by :func:`load_results`.

    Returns:
        The same dictionary, augmented with derived arrays.
    """
    meta = data["metadata"]
    nodes = data["nodes_df"]
    elems = data["elems_df"]

    n_nodes = len(nodes)
    n_elem = len(elems)

    # Pipe section properties
    outer_r = meta["pipe_od_m"] / 2.0  # outer radius [m]
    t_wall = meta["pipe_t_m"]  # wall thickness [m]
    r_mid = outer_r - t_wall / 2.0  # mid-wall radius [m]
    e_mod = meta["pipe_E_Pa"]  # Young's modulus [Pa]

    # A_pipe = 2 * pi * r_mid * t_wall  (thin-walled tube area) [m^2]
    a_pipe = 2.0 * np.pi * r_mid * t_wall

    # I_pipe = pi * r_mid^3 * t_wall  (thin-walled tube moment of inertia) [m^4]
    i_pipe = np.pi * r_mid**3 * t_wall

    # Bending moment at each node
    # Use i-end of element to the right; last node uses j-end of last element
    # M_node[i] = -My_i[i]  (negate: eleForce gives reaction, not internal moment)
    my_i = elems["My_i_Nm"].values  # [N·m]
    my_j = elems["My_j_Nm"].values  # [N·m]
    m_node = np.zeros(n_nodes)
    m_node[:n_elem] = -my_i  # negate reaction sign for internal moment
    m_node[-1] = my_j[-1]  # last node: j-end of last element

    # Axial force at each node (element to the right, last uses last elem)
    n_elem_arr = elems["N_N"].values  # [N]
    n_node = np.zeros(n_nodes)
    n_node[:n_elem] = n_elem_arr
    n_node[-1] = n_elem_arr[-1]

    # sigma_axial = N / A_pipe [Pa]
    sigma_axial = n_node / a_pipe

    # Hoop stress: sigma_hoop = P * R / t  [Pa]
    p_internal = meta.get("P_internal_Pa", 0.0)  # [Pa]
    sigma_hoop = p_internal * outer_r / t_wall if t_wall > 0 else 0.0

    # Determine vertical displacement array
    # For 2D models UZ is all zeros — fall back to UY as "vertical"
    uz = nodes["UZ_m"].values
    uy = nodes["UY_m"].values
    disp_vert = uz if np.any(np.abs(uz) > 1e-15) else uy

    data.update(
        {
            "outer_R": outer_r,
            "t_wall": t_wall,
            "r_mid": r_mid,
            "E": e_mod,
            "A_pipe": a_pipe,
            "I_pipe": i_pipe,
            "M_node": m_node,
            "N_node": n_node,
            "sigma_axial": sigma_axial,
            "sigma_hoop": sigma_hoop,
            "P_internal": p_internal,
            "disp_vert": disp_vert,
        }
    )
    return data


# ── 3D Pipe Surface ──────────────────────────────────────────────────────


def build_pipe_surface(data: dict, deformed: bool = True) -> pv.StructuredGrid:
    """Build a 3D tube surface as a PyVista StructuredGrid.

    At each node, a ring of ``N_CIRC + 1`` points is generated in the
    cross-section plane perpendicular to the local pipe axis.  Adjacent
    rings are connected via the StructuredGrid topology.

    Args:
        data: Data dictionary with node coordinates and pipe geometry.
        deformed: If ``True`` use deformed coordinates, else original.

    Returns:
        PyVista StructuredGrid representing the pipe surface.
    """
    nodes = data["nodes_df"]
    outer_r = data["outer_R"]

    x_orig = nodes["x_m"].values  # [m]
    ux = nodes["UX_m"].values  # [m]
    uy = nodes["UY_m"].values  # [m]
    uz = nodes["UZ_m"].values  # [m]

    n = len(x_orig)

    if deformed:
        px = x_orig + ux  # [m]
        py = uy  # [m]  (original y = 0 for straight pipe)
        pz = uz  # [m]  (original z = 0 for straight pipe)
    else:
        px = x_orig.copy()
        py = np.zeros(n)
        pz = np.zeros(n)

    # ── Local axis directions ─────────────────────────────────────────
    tangents = np.zeros((n, 3))
    tangents[0] = [px[1] - px[0], py[1] - py[0], pz[1] - pz[0]]
    tangents[-1] = [px[-1] - px[-2], py[-1] - py[-2], pz[-1] - pz[-2]]
    for i in range(1, n - 1):
        tangents[i] = [
            px[i + 1] - px[i - 1],
            py[i + 1] - py[i - 1],
            pz[i + 1] - pz[i - 1],
        ]
    norms = np.linalg.norm(tangents, axis=1, keepdims=True)
    norms = np.where(norms < 1e-15, 1.0, norms)
    tangents /= norms

    # ── Build circumferential rings ───────────────────────────────────
    nc = N_CIRC + 1  # +1 to close the tube
    thetas = np.linspace(0, 2 * np.pi, nc)  # [rad]
    cos_t = np.cos(thetas)  # (nc,)
    sin_t = np.sin(thetas)  # (nc,)

    grid_x = np.zeros((n, nc))
    grid_y = np.zeros((n, nc))
    grid_z = np.zeros((n, nc))

    for i in range(n):
        t = tangents[i]

        # Choose initial "up" vector not parallel to tangent
        up = np.array([0.0, 1.0, 0.0])
        if abs(np.dot(t, up)) > 0.99:
            up = np.array([0.0, 0.0, 1.0])

        # Gram-Schmidt: e1 = up − (up·t)*t, normalised
        e1 = up - np.dot(up, t) * t
        e1 /= np.linalg.norm(e1)

        # e2 = t × e1
        e2 = np.cross(t, e1)
        e2 /= np.linalg.norm(e2)

        # point = centre + R * (cos(θ)*e1 + sin(θ)*e2) [m]
        grid_x[i] = px[i] + outer_r * (cos_t * e1[0] + sin_t * e2[0])
        grid_y[i] = py[i] + outer_r * (cos_t * e1[1] + sin_t * e2[1])
        grid_z[i] = pz[i] + outer_r * (cos_t * e1[2] + sin_t * e2[2])

    return pv.StructuredGrid(
        grid_x.reshape(n, nc, 1),
        grid_y.reshape(n, nc, 1),
        grid_z.reshape(n, nc, 1),
    )


# ── Scalar Fields ─────────────────────────────────────────────────────────


def compute_scalar_fields(data: dict) -> dict[str, np.ndarray]:
    """Compute all scalar fields on the pipe surface.

    Returns flattened arrays suitable for assignment to a StructuredGrid's
    ``point_data``.  Theta-varying fields (stress, strain) are computed at
    each circumferential point.  Ring-constant fields (displacement, moment,
    spring force) are repeated around the ring.

    Args:
        data: Data dictionary with derived quantities.

    Returns:
        Dictionary mapping result keys to flattened scalar arrays.
    """
    springs = data["springs_df"]
    nc = N_CIRC + 1

    outer_r = data["outer_R"]
    i_pipe = data["I_pipe"]
    e_mod = data["E"]
    m_node = data["M_node"]
    sigma_axial = data["sigma_axial"]
    sigma_hoop = data["sigma_hoop"]

    thetas = np.linspace(0, 2 * np.pi, nc)  # [rad]
    cos_t = np.cos(thetas)  # (nc,)

    # ── Theta-varying fields (vectorised) ─────────────────────────────
    # sigma_bending(theta) = M * (R * cos(theta)) / I  [Pa]
    sigma_bend = np.outer(m_node, outer_r * cos_t / i_pipe)  # (n, nc)

    # sigma_long = sigma_axial + sigma_bending  [Pa]
    sigma_long = sigma_axial[:, np.newaxis] + sigma_bend  # (n, nc)

    # sigma_vm = sqrt(sL^2 − sL*sH + sH^2)  [Pa] — thin-wall Von Mises
    sigma_vm = np.sqrt(sigma_long**2 - sigma_long * sigma_hoop + sigma_hoop**2)

    # epsilon_long = sigma_long / E  [-]
    epsilon_long = sigma_long / e_mod

    # ── Hoop stress (constant around ring unless pressure varies) ────
    # sigma_hoop = P * R / t  [Pa]  — scalar broadcast to surface
    sigma_hoop_val = data["sigma_hoop"]  # scalar [Pa]
    n_nodes = len(m_node)

    # ── Ring-constant fields ──────────────────────────────────────────
    disp_vert = data["disp_vert"]
    nodes_df = data["nodes_df"]
    disp_ax_nodes = nodes_df["UX_m"].values  # axial displacement [m]

    def ring_constant(vals: np.ndarray) -> np.ndarray:
        """Repeat per-node values around the circumference and flatten."""
        return np.repeat(vals, nc)

    return {
        "sigma_long": sigma_long.ravel(),
        "sigma_vm": sigma_vm.ravel(),
        "sigma_hoop": np.full(n_nodes * nc, sigma_hoop_val),
        "epsilon_long": epsilon_long.ravel(),
        "disp_vert": ring_constant(disp_vert),
        "disp_ax": ring_constant(disp_ax_nodes),
        "moment": ring_constant(m_node),
        "spring_vt": ring_constant(springs["vertical_spring_N"].values),
        "spring_ax": ring_constant(springs["axial_spring_N"].values),
    }


# ── 2D Profile Extraction ────────────────────────────────────────────────


def get_2d_profile(
    data: dict,
    fields: dict[str, np.ndarray],
    result_key: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Extract a 1D profile along the pipe axis for the 2D companion plot.

    For theta-varying fields the extreme-fibre envelope (max absolute value
    between top and bottom) is returned.  For ring-constant fields the nodal
    value is returned directly.

    Args:
        data: Data dictionary.
        fields: Scalar fields from :func:`compute_scalar_fields`.
        result_key: Active result key.

    Returns:
        ``(x_positions [m], values in display units)``.
    """
    x = data["nodes_df"]["x_m"].values
    nc = N_CIRC + 1
    scale = _DISPLAY_SCALE[result_key]
    full = fields[result_key]

    if result_key in ("sigma_long", "sigma_vm", "epsilon_long"):
        # Extract top (theta=0) and bottom (theta=pi) fibres
        vals_top = full[::nc]  # theta = 0
        half = nc // 2
        vals_bot = full[half::nc]  # theta ≈ π
        # Envelope: pick fibre with larger absolute value
        vals = np.where(np.abs(vals_top) >= np.abs(vals_bot), vals_top, vals_bot)
    elif result_key == "sigma_hoop":
        # Hoop stress is constant — just take first of each ring
        vals = full[::nc]
    else:
        vals = full[::nc]

    return x, vals * scale


# ── Yield Detection ───────────────────────────────────────────────────────


def _compute_yield_masks(data: dict) -> tuple[np.ndarray, np.ndarray]:
    """Return boolean masks for yielded vertical and axial springs.

    Uses metadata keys ``vert_yield_per_m_N`` / ``axial_yield_per_m_N``
    (preferred) or falls back to ``k_vt_ul * dy_soil`` /
    ``k_ax_ul * dy_soil_ax``.

    Args:
        data: Data dictionary.

    Returns:
        ``(vt_yielded_mask, ax_yielded_mask)`` — boolean arrays of
        length ``n_nodes``.
    """
    meta = data["metadata"]
    springs = data["springs_df"]
    trib = springs["trib_m"].values  # [m]
    vt_abs = np.abs(springs["vertical_spring_N"].values)  # [N]
    ax_abs = np.abs(springs["axial_spring_N"].values)  # [N]
    n = len(trib)

    # ── Vertical yield ────────────────────────────────────────────────
    vy_per_m = meta.get("vert_yield_per_m_N", None)
    if vy_per_m is None:
        k_vt = meta.get("k_vt_ul", None)
        dy = meta.get("dy_soil", None)
        if k_vt is not None and dy is not None:
            vy_per_m = k_vt * dy  # yield force per unit length [N/m]

    if vy_per_m is not None:
        vt_yield = vy_per_m * trib  # yield force per node [N]
        vt_mask = vt_abs >= 0.99 * vt_yield
    else:
        vt_mask = np.zeros(n, dtype=bool)

    # ── Axial yield ───────────────────────────────────────────────────
    ay_per_m = meta.get("axial_yield_per_m_N", None)
    if ay_per_m is None:
        k_ax = meta.get("k_ax_ul", None)
        dy_ax = meta.get("dy_soil_ax", meta.get("dy_soil", None))
        if k_ax is not None and dy_ax is not None:
            ay_per_m = k_ax * dy_ax  # [N/m]

    if ay_per_m is not None:
        ax_yield = ay_per_m * trib  # [N]
        ax_mask = ax_abs >= 0.99 * ax_yield
    else:
        ax_mask = np.zeros(n, dtype=bool)

    return vt_mask, ax_mask


# ── Terminal Summary ──────────────────────────────────────────────────────


def print_summary(data: dict, fields: dict[str, np.ndarray]) -> None:
    """Print a clean results summary to the terminal.

    Args:
        data: Data dictionary.
        fields: Scalar fields dictionary.
    """
    meta = data["metadata"]
    nodes = data["nodes_df"]
    elems = data["elems_df"]
    n_nodes = len(nodes)
    n_elem = len(elems)

    od_mm = meta["pipe_od_m"] * 1000  # [mm]
    t_mm = meta["pipe_t_m"] * 1000  # [mm]
    e_gpa = meta["pipe_E_Pa"] / 1e9  # [GPa]

    disp_vert = data["disp_vert"]
    max_d_idx = np.argmax(np.abs(disp_vert))
    max_d_mm = disp_vert[max_d_idx] * 1000  # [mm]
    max_d_x = nodes["x_m"].values[max_d_idx]  # [m]

    m_node = data["M_node"]
    max_m_idx = np.argmax(np.abs(m_node))
    max_m = m_node[max_m_idx]  # [N·m]
    max_m_x = nodes["x_m"].values[max_m_idx]

    max_sl = np.max(np.abs(fields["sigma_long"])) / 1e6  # [MPa]
    max_vm = np.max(np.abs(fields["sigma_vm"])) / 1e6  # [MPa]

    vt_mask, ax_mask = _compute_yield_masks(data)

    # Format moment string
    abs_m = abs(max_m)
    if abs_m >= 1e6:
        m_str = f"{max_m / 1e6:.3f} MN·m"
    elif abs_m >= 1e3:
        m_str = f"{max_m / 1e3:.3f} kN·m"
    else:
        m_str = f"{max_m:.3f} N·m"

    print("=" * 60)
    print("OpenSees Pipeline Viewer")
    print("=" * 60)
    print(f"Results loaded from: {data['results_dir']}")
    print(f"Timestamp: {meta.get('timestamp', 'N/A')}")
    print(f"Analysis type: {meta.get('analysis_type', 'N/A')}")
    print(f"Pipe OD: {od_mm:.0f} mm  |  Wall: {t_mm:.1f} mm  |  E: {e_gpa:.0f} GPa")
    print(f"Nodes: {n_nodes}  |  Elements: {n_elem}")
    print()
    print("RESULTS SUMMARY")
    print(f"  Max vertical displacement : {max_d_mm:.1f} mm  at x = {max_d_x:.2f} m")
    print(f"  Max bending moment        : {m_str}  at x = {max_m_x:.2f} m")
    print(f"  Max longitudinal stress   : {max_sl:.2f} MPa")
    print(f"  Max Von Mises stress      : {max_vm:.2f} MPa")
    print(f"  Vertical springs yielded  : {np.sum(vt_mask)} / {n_nodes}")
    print(f"  Axial springs yielded     : {np.sum(ax_mask)} / {n_nodes}")
    print()
    print("KEYBOARD SHORTCUTS")
    print("  1-3  : Stress results (longitudinal, Von Mises, hoop)")
    print("  4    : Strain (longitudinal)")
    print("  5-6  : Displacement (vertical UZ, axial UX)")
    print("  7-9  : Forces (moment, vertical spring, axial spring)")
    print("  T/B  : Top / Bottom view")
    print("  R/L  : Right / Left view")
    print("  F/K  : Front / Back view")
    print("  I    : Isometric view (default)")
    print("  C    : Cycle colormap (jet -> coolwarm -> viridis)")
    print("  S    : Save screenshot")
    print("  Q    : Quit")
    print("=" * 60)


# ── Interactive Viewer ────────────────────────────────────────────────────


def create_viewer(data: dict) -> None:
    """Create and launch the interactive 3D + 2D viewer.

    The PyVista window shows the deformed pipe surface coloured by the
    active scalar field.  A separate Matplotlib window shows the
    corresponding 1D profile plot.  Keyboard shortcuts switch results
    and cycle colormaps.

    Args:
        data: Data dictionary with all loaded and derived data.
    """
    # Build 3D surfaces
    deformed_grid = build_pipe_surface(data, deformed=True)
    undeformed_grid = build_pipe_surface(data, deformed=False)

    # Compute scalar fields and assign to grid
    fields = compute_scalar_fields(data)
    for key, values in fields.items():
        deformed_grid.point_data[key] = values

    # Print terminal summary
    print_summary(data, fields)

    # ── Viewer state ──────────────────────────────────────────────────
    state = {
        "active": RESULT_KEYS[0],
        "cmap_idx": 0,
        "cmaps": ["jet", "coolwarm", "viridis"],
    }

    # ── PyVista plotter ───────────────────────────────────────────────
    plotter = pv.Plotter(title="OpenSees Pipeline Viewer")

    # Static actors: undeformed wireframe
    plotter.add_mesh(
        undeformed_grid,
        style="wireframe",
        color="grey",
        opacity=0.5,
    )

    # Static actors: fault zone (translucent red box)
    meta = data["metadata"]
    if "fault_zone_x_min" in meta and "fault_zone_x_max" in meta:
        x_min = float(meta["fault_zone_x_min"])
        x_max = float(meta["fault_zone_x_max"])
        r = data["outer_R"]
        box = pv.Box(bounds=[x_min, x_max, -3 * r, 3 * r, -3 * r, 3 * r])
        plotter.add_mesh(box, color="red", opacity=0.15, label="Fault Zone")

    # Static actors: yielded spring markers (red spheres)
    vt_mask, _ = _compute_yield_masks(data)
    if np.any(vt_mask):
        nodes_df = data["nodes_df"]
        xs = nodes_df["x_m"].values[vt_mask] + nodes_df["UX_m"].values[vt_mask]
        ys = nodes_df["UY_m"].values[vt_mask] - data["outer_R"]
        zs = nodes_df["UZ_m"].values[vt_mask]
        sphere_pts = np.column_stack([xs, ys, zs])
        plotter.add_mesh(
            pv.PolyData(sphere_pts),
            color="red",
            point_size=10,
            render_points_as_spheres=True,
            label="Yielded Springs",
        )

    # Dynamic actor: deformed pipe (use name= for replacement)
    def _add_deformed() -> None:
        """Add (or replace) the deformed pipe mesh with current settings."""
        info = RESULT_INFO[state["active"]]
        cmap = state["cmaps"][state["cmap_idx"]]
        plotter.add_mesh(
            deformed_grid,
            scalars=state["active"],
            cmap=cmap,
            name="deformed_pipe",
            scalar_bar_args={
                "title": f"{info['label']} [{info['unit']}]",
                "shadow": True,
            },
        )

    _add_deformed()

    # ── On-screen menu overlay ────────────────────────────────────────
    _MENU_BODY = (
        "\n"
        "[STRESS]           [STRAIN]        [DISPLACEMENT]    [FORCES]\n"
        "1 Longitudinal     4 Long. Strain   5 Vertical UZ     7 Moment\n"
        "2 Von Mises                         6 Axial UX        8 Vert. Spring\n"
        "3 Hoop                                                9 Axial Spring\n"
        "\n"
        "[VIEWS]\n"
        "T Top  B Bottom  R Right  L Left  F Front  K Back  I Iso\n"
        "\n"
        "C Cycle colormap   S Save screenshot   Q Quit"
    )

    # Mutable reference for the header text actor so we can remove/re-add it
    state["menu_actor"] = None

    def _update_menu() -> None:
        """Update the on-screen result header + static menu overlay."""
        if state["menu_actor"] is not None:
            plotter.remove_actor(state["menu_actor"])
        info = RESULT_INFO[state["active"]]
        header = f"RESULT: {info['label']} ({info['unit']})"
        full_text = header + _MENU_BODY
        state["menu_actor"] = plotter.add_text(
            full_text,
            position="upper_right",
            font_size=8,
            font="courier",
            color="white",
            name="menu_overlay",
        )

    _update_menu()

    # ── Default camera ────────────────────────────────────────────────
    plotter.view_isometric()

    # ── Matplotlib 2D companion ───────────────────────────────────────
    plt.ion()
    fig_2d, ax_2d = plt.subplots(figsize=(6, 4))
    fig_2d.canvas.manager.set_window_title("Pipeline Profile")

    def _update_2d() -> None:
        """Redraw the 2D companion plot for the active result."""
        ax_2d.clear()
        rk = state["active"]
        x, vals = get_2d_profile(data, fields, rk)
        info = RESULT_INFO[rk]
        ax_2d.plot(x, vals, "b-", linewidth=1.5)
        ax_2d.set_xlabel("Distance along pipe [m]")
        ax_2d.set_ylabel(_DISPLAY_YLABEL[rk])
        ax_2d.set_title(info["label"])
        ax_2d.grid(True, alpha=0.3)
        fig_2d.tight_layout()
        fig_2d.canvas.draw_idle()
        fig_2d.canvas.flush_events()

    _update_2d()

    # ── Key bindings: result switching (1–9) ──────────────────────────

    def _switch(key: str) -> None:
        idx = int(key) - 1
        if 0 <= idx < len(RESULT_KEYS):
            result_key = RESULT_KEYS[idx]
            # Guard: do not switch to empty or all-NaN fields
            values = deformed_grid.point_data[result_key]
            if np.all(values == 0) or np.all(np.isnan(values)):
                print(
                    f"  WARNING: {result_key} has no data"
                    " — staying on current result"
                )
                return
            state["active"] = result_key
            _add_deformed()
            _update_menu()
            plotter.render()
            _update_2d()
            print(f"  → {RESULT_INFO[state['active']]['label']}")

    for k in "123456789":
        plotter.add_key_event(k, lambda key=k: _switch(key))

    # ── Key bindings: camera views ────────────────────────────────────
    plotter.add_key_event("t", lambda: plotter.view_xy())
    plotter.add_key_event("b", lambda: plotter.view_xy(negative=True))
    plotter.add_key_event("r", lambda: plotter.view_yz())
    plotter.add_key_event("l", lambda: plotter.view_yz(negative=True))
    plotter.add_key_event("f", lambda: plotter.view_xz())
    plotter.add_key_event("k", lambda: plotter.view_xz(negative=True))
    plotter.add_key_event("i", lambda: plotter.view_isometric())

    # ── Key bindings: colormap cycling ────────────────────────────────

    def _cycle_cmap() -> None:
        state["cmap_idx"] = (state["cmap_idx"] + 1) % len(state["cmaps"])
        cmap = state["cmaps"][state["cmap_idx"]]
        print(f"  → Colormap: {cmap}")
        _add_deformed()
        plotter.render()

    plotter.add_key_event("c", _cycle_cmap)

    # ── Key bindings: screenshot ──────────────────────────────────────

    def _screenshot() -> None:
        rk = state["active"]
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        fname = f"viewer_screenshot_{rk}_{ts}.png"
        fpath = Path(data["results_dir"]) / fname
        plotter.screenshot(str(fpath))
        print(f"  Screenshot saved: {fpath}")

    plotter.add_key_event("s", _screenshot)

    # ── Show ──────────────────────────────────────────────────────────
    plotter.show()
    plt.close(fig_2d)


# ── Entry Point ───────────────────────────────────────────────────────────


def main() -> None:
    """Entry point for the standalone pipeline viewer."""
    # Handle drag-and-drop path input
    if len(sys.argv) > 1:
        results_dir = sys.argv[1].strip()
    else:
        results_dir = input(
            "Drag and drop your results folder here, then press Enter: "
        ).strip()

    # Strip quotes that Mac sometimes adds around paths with spaces
    results_dir = results_dir.strip("'\"")

    if not os.path.isdir(results_dir):
        print(f"ERROR: Directory not found: {results_dir}")
        sys.exit(1)

    data = load_results(results_dir)
    data = compute_derived(data)
    create_viewer(data)


if __name__ == "__main__":
    main()
