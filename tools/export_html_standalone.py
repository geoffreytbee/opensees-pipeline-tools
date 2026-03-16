#!/usr/bin/env python3
"""Standalone HTML exporter for OpenSees pipeline results.

Reads a results folder produced by ``opensees_pipeline.postprocess.write_results``
and writes a single self-contained interactive HTML file using Plotly for 3D
rendering.  No server, internet, or Python installation is required to view
the output — open it in any modern browser.

Usage::

    python tools/export_html_standalone.py results/benchmark_01
    python tools/export_html_standalone.py          # drag-and-drop prompt

Dependencies: numpy, pandas, plotly.
"""

from __future__ import annotations

import json
import os
import sys
from datetime import datetime
from typing import Any

import numpy as np
import pandas as pd

try:
    import plotly.graph_objects as go
    import plotly.offline
except ImportError:
    print("ERROR: plotly is not installed. Run: pip install plotly")
    sys.exit(1)

# ── Constants ─────────────────────────────────────────────────────────────

N_CIRC = 16  # circumferential discretisation points per ring

DEFORM_SCALES = [0, 1, 5, 10, 25, 50, 100]

SCALAR_FIELDS: dict[str, dict[str, str]] = {
    "sigma_long": {"label": "Longitudinal Stress", "unit": "Pa"},
    "sigma_vm": {"label": "Von Mises Stress", "unit": "Pa"},
    "sigma_hoop": {"label": "Hoop Stress", "unit": "Pa"},
    "epsilon_long": {"label": "Longitudinal Strain", "unit": "-"},
    "disp_uz": {"label": "Vertical Displacement", "unit": "m"},
    "disp_ux": {"label": "Axial Displacement", "unit": "m"},
    "moment": {"label": "Bending Moment", "unit": "N·m"},
    "spring_vt": {"label": "Vertical Spring Force", "unit": "N"},
    "spring_ax": {"label": "Axial Spring Force", "unit": "N"},
}


# ── Data Loading ──────────────────────────────────────────────────────────


def load_results(results_dir: str) -> dict[str, Any]:
    """Load all results files from *results_dir*.

    Args:
        results_dir: Path to directory containing the four standard files.

    Returns:
        Dictionary with DataFrames and metadata.
    """
    rd = results_dir
    nodes_df = pd.read_csv(os.path.join(rd, "pipe_nodes.csv"))
    elems_df = pd.read_csv(os.path.join(rd, "pipe_elements.csv"))
    springs_df = pd.read_csv(os.path.join(rd, "spring_forces.csv"))
    with open(os.path.join(rd, "model_metadata.json")) as fh:
        metadata = json.load(fh)
    return {
        "nodes_df": nodes_df,
        "elems_df": elems_df,
        "springs_df": springs_df,
        "metadata": metadata,
    }


# ── Derived Quantities ────────────────────────────────────────────────────


def compute_derived(data: dict[str, Any]) -> dict[str, Any]:
    """Compute derived section properties and per-node stress arrays.

    Args:
        data: Dictionary from :func:`load_results`.

    Returns:
        Updated dictionary with derived arrays.
    """
    meta = data["metadata"]
    nodes = data["nodes_df"]
    elems = data["elems_df"]
    n_nodes = len(nodes)
    n_elem = len(elems)

    outer_r = meta["pipe_od_m"] / 2.0  # [m]
    t_wall = meta["pipe_t_m"]  # [m]
    r_mid = outer_r - t_wall / 2.0  # [m]
    e_mod = meta["pipe_E_Pa"]  # [Pa]

    # A = 2 * pi * r_mid * t  [m^2]
    a_pipe = 2.0 * np.pi * r_mid * t_wall
    # I = pi * r_mid^3 * t  [m^4]
    i_pipe = np.pi * r_mid**3 * t_wall

    # Bending moment at each node: -My_i of element to the right
    my_i = elems["My_i_Nm"].values  # [N·m]
    my_j = elems["My_j_Nm"].values  # [N·m]
    m_node = np.zeros(n_nodes)
    m_node[:n_elem] = -my_i
    m_node[-1] = my_j[-1]

    # Axial force at each node
    n_arr = elems["N_N"].values  # [N]
    n_node = np.zeros(n_nodes)
    n_node[:n_elem] = n_arr
    n_node[-1] = n_arr[-1]

    # sigma_axial = N / A  [Pa]
    sigma_axial = n_node / a_pipe

    # sigma_hoop = P * R / t  [Pa]
    p_int = meta.get("P_internal_Pa", 0.0)  # [Pa]
    sigma_hoop = p_int * outer_r / t_wall if t_wall > 0 else 0.0

    # Vertical displacement — use UZ, fall back to UY for 2D models
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
            "sigma_hoop_scalar": sigma_hoop,
            "disp_vert": disp_vert,
        }
    )
    return data


# ── Scalar Fields ─────────────────────────────────────────────────────────


def compute_scalar_fields(data: dict[str, Any]) -> dict[str, np.ndarray]:
    """Compute all scalar fields at every surface point.

    Each field has shape ``(n_nodes * (N_CIRC + 1),)``.

    Args:
        data: Dictionary with derived quantities.

    Returns:
        Mapping from field key to flat float array.
    """
    nodes = data["nodes_df"]
    springs = data["springs_df"]
    nc = N_CIRC + 1
    n_nodes = len(nodes)

    outer_r = data["outer_R"]
    i_pipe = data["I_pipe"]
    e_mod = data["E"]
    m_node = data["M_node"]
    sigma_axial = data["sigma_axial"]
    sigma_hoop = data["sigma_hoop_scalar"]

    thetas = np.linspace(0, 2 * np.pi, nc)  # [rad]
    cos_t = np.cos(thetas)

    # sigma_bending(theta) = M * R * cos(theta) / I  [Pa]
    sigma_bend = np.outer(m_node, outer_r * cos_t / i_pipe)

    # sigma_long = sigma_axial + sigma_bending  [Pa]
    sigma_long = sigma_axial[:, np.newaxis] + sigma_bend

    # sigma_vm = sqrt(sL^2 - sL*sH + sH^2)  [Pa]
    sigma_vm = np.sqrt(sigma_long**2 - sigma_long * sigma_hoop + sigma_hoop**2)

    # epsilon_long = sigma_long / E  [-]
    epsilon_long = sigma_long / e_mod

    def ring_constant(vals: np.ndarray) -> np.ndarray:
        """Repeat per-node values around circumference."""
        return np.repeat(vals, nc)

    return {
        "sigma_long": sigma_long.ravel(),
        "sigma_vm": sigma_vm.ravel(),
        "sigma_hoop": np.full(n_nodes * nc, sigma_hoop),
        "epsilon_long": epsilon_long.ravel(),
        "disp_uz": ring_constant(data["disp_vert"]),
        "disp_ux": ring_constant(nodes["UX_m"].values),
        "moment": ring_constant(m_node),
        "spring_vt": ring_constant(springs["vertical_spring_N"].values),
        "spring_ax": ring_constant(springs["axial_spring_N"].values),
    }


# ── Tube Geometry ─────────────────────────────────────────────────────────


def generate_tube_mesh(
    x_nodes: np.ndarray,
    disp_ux: np.ndarray,
    disp_uy: np.ndarray,
    disp_uz: np.ndarray,
    outer_r: float,
    deform_scale: float,
    geom_scale: float,
    radial_scale: float,
    n_circ: int = N_CIRC,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Generate tube surface as triangle mesh arrays for Plotly Mesh3d.

    Coordinates are normalised: axial by *geom_scale*, radial by
    *geom_scale * radial_scale*.

    Args:
        x_nodes: Original node X-positions [m].
        disp_ux: Nodal X-displacements [m].
        disp_uy: Nodal Y-displacements [m].
        disp_uz: Nodal Z-displacements [m].
        outer_r: Pipe outer radius [m].
        deform_scale: Deformation magnification factor [-].
        geom_scale: Axial normalisation factor [1/m].
        radial_scale: Radial exaggeration factor [-].
        n_circ: Circumferential point count.

    Returns:
        ``(x, y, z, i_tri, j_tri, k_tri)`` arrays for Plotly Mesh3d.
    """
    n = len(x_nodes)
    nc = n_circ + 1  # +1 to close tube

    # Deformed centre-line [m]
    px = x_nodes + disp_ux * deform_scale
    py = disp_uy * deform_scale
    pz = disp_uz * deform_scale

    # Tangent vectors
    tangents = np.zeros((n, 3))
    tangents[0] = [px[1] - px[0], py[1] - py[0], pz[1] - pz[0]]
    tangents[-1] = [px[-1] - px[-2], py[-1] - py[-2], pz[-1] - pz[-2]]
    for idx in range(1, n - 1):
        tangents[idx] = [
            px[idx + 1] - px[idx - 1],
            py[idx + 1] - py[idx - 1],
            pz[idx + 1] - pz[idx - 1],
        ]
    norms = np.linalg.norm(tangents, axis=1, keepdims=True)
    norms = np.where(norms < 1e-15, 1.0, norms)
    tangents /= norms

    thetas = np.linspace(0, 2 * np.pi, nc)
    cos_t = np.cos(thetas)
    sin_t = np.sin(thetas)

    # Build surface points [m] then normalise
    pts = np.zeros((n * nc, 3))
    for idx in range(n):
        t = tangents[idx]
        up = np.array([0.0, 1.0, 0.0])
        if abs(np.dot(t, up)) > 0.99:
            up = np.array([0.0, 0.0, 1.0])
        e1 = up - np.dot(up, t) * t
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(t, e1)
        e2 /= np.linalg.norm(e2)

        base = idx * nc
        pts[base : base + nc, 0] = px[idx] + outer_r * (cos_t * e1[0] + sin_t * e2[0])
        pts[base : base + nc, 1] = py[idx] + outer_r * (cos_t * e1[1] + sin_t * e2[1])
        pts[base : base + nc, 2] = pz[idx] + outer_r * (cos_t * e1[2] + sin_t * e2[2])

    # Normalise coordinates for browser rendering
    pts[:, 0] *= geom_scale  # axial
    pts[:, 1] *= geom_scale * radial_scale  # radial y
    pts[:, 2] *= geom_scale * radial_scale  # radial z

    # Build triangle connectivity (two triangles per quad)
    i_list: list[int] = []
    j_list: list[int] = []
    k_list: list[int] = []
    for row in range(n - 1):
        for col in range(n_circ):
            col_next = (col + 1) % nc
            p0 = row * nc + col
            p1 = row * nc + col_next
            p2 = (row + 1) * nc + col_next
            p3 = (row + 1) * nc + col
            # Triangle 1: p0-p1-p2
            i_list.append(p0)
            j_list.append(p1)
            k_list.append(p2)
            # Triangle 2: p0-p2-p3
            i_list.append(p0)
            j_list.append(p2)
            k_list.append(p3)

    return (
        pts[:, 0],
        pts[:, 1],
        pts[:, 2],
        np.array(i_list),
        np.array(j_list),
        np.array(k_list),
    )


# ── Yield Detection ──────────────────────────────────────────────────────


def _yielded_nodes(data: dict[str, Any]) -> list[int]:
    """Return indices of nodes with yielded vertical springs.

    Args:
        data: Data dictionary.

    Returns:
        List of 0-based node indices.
    """
    meta = data["metadata"]
    springs = data["springs_df"]
    trib = springs["trib_m"].values
    vt_abs = np.abs(springs["vertical_spring_N"].values)

    vy = meta.get("vert_yield_per_m_N", None)
    if vy is None:
        k = meta.get("k_vt_ul", None)
        dy = meta.get("dy_soil", None)
        if k is not None and dy is not None:
            vy = k * dy
    if vy is None:
        return []

    vt_yield = vy * trib
    return [int(i) for i in np.where(vt_abs >= 0.99 * vt_yield)[0]]


# ── Field Statistics ──────────────────────────────────────────────────────


def compute_field_stats(
    fields: dict[str, np.ndarray],
    x_nodes: np.ndarray,
    n_nodes: int,
) -> dict[str, dict[str, Any]]:
    """Compute per-field min/max values and their x-locations.

    Args:
        fields: Scalar fields dictionary.
        x_nodes: Node x-positions [m].
        n_nodes: Number of pipe nodes.

    Returns:
        Dictionary mapping field keys to stat dicts.
    """
    nc = N_CIRC + 1
    stats: dict[str, dict[str, Any]] = {}
    for key, info in SCALAR_FIELDS.items():
        vals = fields[key]
        available = not (np.all(vals == 0) or np.all(np.isnan(vals)))

        # Per-node envelope (max absolute across circumference)
        node_vals = np.array(
            [
                vals[i * nc : (i + 1) * nc][
                    np.argmax(np.abs(vals[i * nc : (i + 1) * nc]))
                ]
                for i in range(n_nodes)
            ]
        )
        max_idx = int(np.argmax(node_vals))
        min_idx = int(np.argmin(node_vals))

        stats[key] = {
            "min": float(node_vals[min_idx]),
            "max": float(node_vals[max_idx]),
            "min_x": float(x_nodes[min_idx]),
            "max_x": float(x_nodes[max_idx]),
            "unit": info["unit"],
            "label": info["label"],
            "available": available,
        }
    return stats


# ── Plotly Figure Builder ─────────────────────────────────────────────────


def build_figure(
    data: dict[str, Any],
    fields: dict[str, np.ndarray],
    field_stats: dict[str, dict[str, Any]],
    default_scale: int = 10,
) -> go.Figure:
    """Build a Plotly figure with all result traces and dropdown controls.

    Each available scalar field gets its own Mesh3d trace (initially hidden
    except for the first).  A dropdown selector toggles visibility.

    Args:
        data: Data dictionary with derived quantities.
        fields: Scalar fields dictionary.
        field_stats: Per-field statistics.
        default_scale: Default deformation scale factor.

    Returns:
        Plotly Figure ready for export.
    """
    nodes = data["nodes_df"]
    meta = data["metadata"]
    x_nodes = nodes["x_m"].values
    outer_r = data["outer_R"]
    n_nodes = len(x_nodes)
    nc = N_CIRC + 1

    pipe_length = float(x_nodes[-1] - x_nodes[0])  # [m]
    geom_scale = 1.0 / pipe_length if pipe_length > 0 else 1.0
    radial_scale = 20.0

    ux = nodes["UX_m"].values
    uy = nodes["UY_m"].values
    uz = nodes["UZ_m"].values

    # Deformed mesh geometry
    xm, ym, zm, ii, jj, kk = generate_tube_mesh(
        x_nodes,
        ux,
        uy,
        uz,
        outer_r,
        float(default_scale),
        geom_scale,
        radial_scale,
    )

    # Undeformed wireframe geometry
    xu, yu, zu, _, _, _ = generate_tube_mesh(
        x_nodes,
        np.zeros_like(ux),
        np.zeros_like(uy),
        np.zeros_like(uz),
        outer_r,
        0.0,
        geom_scale,
        radial_scale,
    )

    # Collect available field keys
    avail_keys = [k for k, s in field_stats.items() if s["available"]]
    n_avail = len(avail_keys)

    fig = go.Figure()

    # ── Add one Mesh3d trace per available field ──────────────────────
    for trace_idx, fkey in enumerate(avail_keys):
        info = field_stats[fkey]
        # Per-vertex intensity (one value per surface point)
        intensity = fields[fkey].astype(np.float64)

        fig.add_trace(
            go.Mesh3d(
                x=xm,
                y=ym,
                z=zm,
                i=ii,
                j=jj,
                k=kk,
                intensity=intensity,
                intensitymode="vertex",
                colorscale="Jet",
                colorbar=dict(
                    title=dict(
                        text=f"{info['label']}<br>({info['unit']})",
                        side="right",
                    ),
                    len=0.75,
                ),
                cmin=info["min"],
                cmax=info["max"],
                flatshading=True,
                lighting=dict(ambient=1.0, diffuse=0, specular=0),
                visible=(trace_idx == 0),
                name=info["label"],
                hovertemplate=(
                    f"{info['label']}: %{{intensity:.3e}} {info['unit']}"
                    "<extra></extra>"
                ),
            )
        )

    # ── Undeformed wireframe (always visible) ─────────────────────────
    # Sample lines along axial and circumferential directions
    wire_x: list[float | None] = []
    wire_y: list[float | None] = []
    wire_z: list[float | None] = []

    # Axial lines (every 4th circumferential point)
    for j in range(0, nc, max(1, nc // 4)):
        for row in range(n_nodes):
            idx = row * nc + j
            wire_x.append(float(xu[idx]))
            wire_y.append(float(yu[idx]))
            wire_z.append(float(zu[idx]))
        wire_x.append(None)
        wire_y.append(None)
        wire_z.append(None)

    # Circumferential rings (every 5th node)
    for row in range(0, n_nodes, max(1, n_nodes // 10)):
        for j in range(nc):
            idx = row * nc + j
            wire_x.append(float(xu[idx]))
            wire_y.append(float(yu[idx]))
            wire_z.append(float(zu[idx]))
        wire_x.append(None)
        wire_y.append(None)
        wire_z.append(None)

    fig.add_trace(
        go.Scatter3d(
            x=wire_x,
            y=wire_y,
            z=wire_z,
            mode="lines",
            line=dict(color="grey", width=1),
            name="Undeformed",
            hoverinfo="skip",
            visible=True,
        )
    )

    # ── Fault zone box ────────────────────────────────────────────────
    fz_min = meta.get("fault_zone_x_min")
    fz_max = meta.get("fault_zone_x_max")
    if fz_min is not None and fz_max is not None:
        fz_min_s = fz_min * geom_scale
        fz_max_s = fz_max * geom_scale
        r_s = outer_r * geom_scale * radial_scale * 1.5
        # 8 corners of the box
        bx = [
            fz_min_s,
            fz_max_s,
            fz_max_s,
            fz_min_s,
            fz_min_s,
            fz_max_s,
            fz_max_s,
            fz_min_s,
        ]
        by = [-r_s, -r_s, r_s, r_s, -r_s, -r_s, r_s, r_s]
        bz = [-r_s, -r_s, -r_s, -r_s, r_s, r_s, r_s, r_s]
        fig.add_trace(
            go.Mesh3d(
                x=bx,
                y=by,
                z=bz,
                i=[0, 0, 0, 0, 4, 4, 2, 2, 0, 0, 1, 1],
                j=[1, 2, 4, 5, 5, 6, 3, 7, 1, 4, 2, 6],
                k=[2, 3, 5, 6, 6, 7, 7, 6, 5, 7, 6, 5],
                color="red",
                opacity=0.12,
                name="Fault Zone",
                hoverinfo="skip",
                visible=True,
                flatshading=True,
            )
        )

    # ── Yielded spring markers ────────────────────────────────────────
    yielded = _yielded_nodes(data)
    if yielded:
        yx = (x_nodes[yielded] + ux[yielded] * default_scale) * geom_scale
        yy = (uy[yielded] * default_scale) * geom_scale * radial_scale
        yz_pts = (uz[yielded] * default_scale - outer_r) * geom_scale * radial_scale
        fig.add_trace(
            go.Scatter3d(
                x=yx,
                y=yy,
                z=yz_pts,
                mode="markers",
                marker=dict(size=4, color="red", symbol="circle"),
                name="Yielded Springs",
                hoverinfo="text",
                text=[f"Node {i}" for i in yielded],
                visible=True,
            )
        )

    # ── Dropdown for result switching ─────────────────────────────────
    # Build visibility arrays for each dropdown option.
    # Traces: [mesh_0, mesh_1, ..., mesh_{n-1}, wireframe, (fault), (yielded)]
    n_extra = 1  # wireframe
    if fz_min is not None:
        n_extra += 1
    if yielded:
        n_extra += 1

    buttons = []
    for btn_idx, fkey in enumerate(avail_keys):
        info = field_stats[fkey]
        vis = [False] * n_avail + [True] * n_extra
        vis[btn_idx] = True
        buttons.append(
            dict(
                label=info["label"],
                method="update",
                args=[{"visible": vis}],
            )
        )

    # ── Layout ────────────────────────────────────────────────────────
    fig.update_layout(
        title=dict(
            text="OpenSees Pipeline Report",
            font=dict(size=16, color="#e0e0e0"),
        ),
        paper_bgcolor="#1a1a2e",
        plot_bgcolor="#1a1a2e",
        font=dict(color="#e0e0e0"),
        scene=dict(
            xaxis=dict(title="Axial", showgrid=True, gridcolor="#333", zeroline=False),
            yaxis=dict(title="Y", showgrid=True, gridcolor="#333", zeroline=False),
            zaxis=dict(title="Z", showgrid=True, gridcolor="#333", zeroline=False),
            bgcolor="#16213e",
            aspectmode="data",
            camera=dict(
                eye=dict(x=1.5, y=-1.5, z=0.8),
                up=dict(x=0, y=0, z=1),
            ),
        ),
        updatemenus=[
            dict(
                buttons=buttons,
                direction="down",
                showactive=True,
                x=0.0,
                xanchor="left",
                y=1.12,
                yanchor="top",
                bgcolor="#0f3460",
                font=dict(color="#e0e0e0"),
            ),
            # Camera view presets
            dict(
                type="buttons",
                buttons=[
                    dict(
                        label="Iso",
                        method="relayout",
                        args=[
                            {
                                "scene.camera": dict(
                                    eye=dict(x=1.5, y=-1.5, z=0.8),
                                    up=dict(x=0, y=0, z=1),
                                )
                            }
                        ],
                    ),
                    dict(
                        label="Top",
                        method="relayout",
                        args=[
                            {
                                "scene.camera": dict(
                                    eye=dict(x=0, y=0, z=3), up=dict(x=0, y=1, z=0)
                                )
                            }
                        ],
                    ),
                    dict(
                        label="Front",
                        method="relayout",
                        args=[
                            {
                                "scene.camera": dict(
                                    eye=dict(x=3, y=0, z=0), up=dict(x=0, y=0, z=1)
                                )
                            }
                        ],
                    ),
                    dict(
                        label="Right",
                        method="relayout",
                        args=[
                            {
                                "scene.camera": dict(
                                    eye=dict(x=0, y=-3, z=0), up=dict(x=0, y=0, z=1)
                                )
                            }
                        ],
                    ),
                ],
                direction="left",
                showactive=False,
                x=0.35,
                xanchor="left",
                y=1.12,
                yanchor="top",
                bgcolor="#0f3460",
                font=dict(color="#e0e0e0"),
            ),
        ],
        legend=dict(
            x=1.0,
            y=0.5,
            bgcolor="rgba(22,33,62,0.8)",
            font=dict(size=11),
        ),
        margin=dict(l=0, r=0, t=80, b=0),
        height=700,
    )

    return fig


# ── Summary Table HTML ────────────────────────────────────────────────────


def build_summary_html(
    field_stats: dict[str, dict[str, Any]],
    meta: dict[str, Any],
) -> str:
    """Build an HTML summary table string to append below the Plotly figure.

    Args:
        field_stats: Per-field statistics.
        meta: Model metadata dictionary.

    Returns:
        HTML string with styled summary table.
    """
    od_mm = meta["pipe_od_mm"]
    t_mm = meta["pipe_t_mm"]
    length_m = meta["pipe_length_m"]
    n_nodes = meta["n_nodes"]
    n_elem = meta["n_elements"]

    rows = ""
    for key, st in field_stats.items():
        if not st["available"]:
            continue
        rows += (
            f"<tr>"
            f"<td>{st['label']}</td>"
            f"<td>{st['max']:.4e}</td>"
            f"<td>{st['max_x']:.2f} m</td>"
            f"<td>{st['min']:.4e}</td>"
            f"<td>{st['min_x']:.2f} m</td>"
            f"<td>{st['unit']}</td>"
            f"</tr>\n"
        )

    return f"""
<div style="background:#16213e;color:#e0e0e0;padding:12px 20px;
            font-family:'Segoe UI',Arial,sans-serif;font-size:13px">
  <p style="margin:4px 0"><b>Pipe:</b> OD={od_mm:.0f} mm,
     t={t_mm:.1f} mm, L={length_m:.1f} m &nbsp;|&nbsp;
     <b>Nodes:</b> {n_nodes} &nbsp;|&nbsp;
     <b>Elements:</b> {n_elem}</p>
  <table style="width:100%;border-collapse:collapse;margin-top:8px">
    <thead>
      <tr style="background:#0f3460">
        <th style="padding:5px 10px;text-align:left">Result</th>
        <th style="padding:5px 10px;text-align:left">Max Value</th>
        <th style="padding:5px 10px;text-align:left">Max @ x</th>
        <th style="padding:5px 10px;text-align:left">Min Value</th>
        <th style="padding:5px 10px;text-align:left">Min @ x</th>
        <th style="padding:5px 10px;text-align:left">Unit</th>
      </tr>
    </thead>
    <tbody>
      {rows}
    </tbody>
  </table>
</div>
"""


# ── Entry Point ───────────────────────────────────────────────────────────


def main() -> None:
    """Entry point for the standalone HTML exporter."""
    if len(sys.argv) > 1:
        results_dir = sys.argv[1].strip().strip("'\"")
    else:
        results_dir = input("Drag results folder here, then press Enter: ").strip()
        results_dir = results_dir.strip("'\"")

    if not os.path.isdir(results_dir):
        print(f"ERROR: Directory not found: {results_dir}")
        sys.exit(1)

    print("=" * 60)
    print("OpenSees Pipeline — Standalone HTML Exporter")
    print("=" * 60)
    print(f"Loading results from: {results_dir}/")

    data = load_results(results_dir)
    data = compute_derived(data)
    nodes = data["nodes_df"]
    meta = data["metadata"]
    x_nodes = nodes["x_m"].values
    n_nodes = len(x_nodes)
    n_elem = len(data["elems_df"])
    od_mm = meta["pipe_od_m"] * 1000  # [mm]
    t_mm = meta["pipe_t_m"] * 1000  # [mm]

    print(f"  Nodes: {n_nodes} | Elements: {n_elem}")
    print(f"  Pipe: OD={od_mm:.0f}mm, t={t_mm:.1f}mm")

    print("\nComputing scalar fields...")
    fields = compute_scalar_fields(data)
    field_stats = compute_field_stats(fields, x_nodes, n_nodes)

    for key, info in SCALAR_FIELDS.items():
        st = field_stats[key]
        if st["available"]:
            print(
                f"  OK  {info['label']:30s} "
                f"min={st['min']:.3e}  max={st['max']:.3e}"
            )
        else:
            print(f"  --  {info['label']:30s} no data")

    print("\nBuilding Plotly figure...")
    fig = build_figure(data, fields, field_stats)

    summary_html = build_summary_html(
        field_stats,
        {
            "pipe_od_mm": od_mm,
            "pipe_t_mm": t_mm,
            "pipe_length_m": float(x_nodes[-1] - x_nodes[0]),
            "n_nodes": n_nodes,
            "n_elements": n_elem,
        },
    )

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    fname = f"pipeline_report_{timestamp}.html"
    output_path = os.path.join(results_dir, fname)

    print("\nWriting output...")
    plotly.offline.plot(
        fig,
        filename=output_path,
        auto_open=False,
        include_plotlyjs=True,
    )

    # Append summary table to the HTML file
    with open(output_path, "a") as fh:
        fh.write(summary_html)

    file_mb = os.path.getsize(output_path) / 1e6
    print(f"  {fname} ({file_mb:.1f} MB)")

    if file_mb > 50:
        print("  WARNING: File is large — may be slow to open in browser")

    print()
    print("Open in any modern browser — no Python required.")
    print("Share by email, Dropbox, or Google Drive.")
    print("=" * 60)


if __name__ == "__main__":
    main()
