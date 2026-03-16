#!/usr/bin/env python3
"""OpenSees Pipeline Viewer — Web UI powered by Panel + PyVista.

A browser-based interactive post-processor for pipeline analysis results.
Opens a browser tab with a professional engineering UI featuring a 3D PyVista
render embedded via Panel's VTK pane, toolbar controls for result switching,
colormap selection, deformation scaling, and standard camera views.

Dependencies: panel, pyvista, numpy, pandas
Does **not** require openseespy.

Usage::

    python tools/opensees_pipeline_viewer_panel.py [results_dir]

If *results_dir* is not provided, the script prompts for drag-and-drop input.
"""

from __future__ import annotations

import os
import sys
from datetime import datetime
from pathlib import Path

import nest_asyncio
import numpy as np
import pandas as pd

nest_asyncio.apply()

# ── Dependency check (panel / pyvista may not be installed) ───────────────
try:
    import panel as pn
    import pyvista as pv
except ImportError as exc:
    print(f"ERROR: {exc}. Run: pip install panel pyvista trame trame-vtk")
    sys.exit(1)

# Import data loading and computation from the existing desktop viewer.
# This avoids duplicating ~300 lines of validated computation code.
_tools_dir = os.path.dirname(os.path.abspath(__file__))
if _tools_dir not in sys.path:
    sys.path.insert(0, _tools_dir)

from opensees_pipeline_viewer import (  # noqa: E402
    N_CIRC,
    _compute_yield_masks,
    build_pipe_surface,
    compute_derived,
    compute_scalar_fields,
    load_results,
    print_summary,
)

# ── Result Definitions ────────────────────────────────────────────────────

RESULTS: dict[str, dict[str, str]] = {
    # STRESS
    "Longitudinal Stress": {"key": "sigma_long", "unit": "Pa", "group": "Stress"},
    "Von Mises Stress": {"key": "sigma_vm", "unit": "Pa", "group": "Stress"},
    "Hoop Stress": {"key": "sigma_hoop", "unit": "Pa", "group": "Stress"},
    # STRAIN
    "Longitudinal Strain": {"key": "epsilon_long", "unit": "-", "group": "Strain"},
    # DISPLACEMENT
    "Vertical UZ": {"key": "disp_vert", "unit": "m", "group": "Displacement"},
    "Axial UX": {"key": "disp_ax", "unit": "m", "group": "Displacement"},
    # FORCES
    "Bending Moment": {"key": "moment", "unit": "N·m", "group": "Forces"},
    "Vertical Spring Force": {"key": "spring_vt", "unit": "N", "group": "Forces"},
    "Axial Spring Force": {"key": "spring_ax", "unit": "N", "group": "Forces"},
}


# ── Scaled Surface Builder ────────────────────────────────────────────────


def _build_scaled_surface(data: dict, scale: float) -> pv.StructuredGrid:
    """Build deformed tube surface with scaled displacements.

    Temporarily scales the displacement columns in the data DataFrame,
    calls ``build_pipe_surface``, then restores the originals.

    Args:
        data: Data dictionary from :func:`compute_derived`.
        scale: Deformation magnification factor [-].

    Returns:
        PyVista StructuredGrid of the deformed tube surface.
    """
    nodes = data["nodes_df"]
    orig_ux = nodes["UX_m"].values.copy()
    orig_uy = nodes["UY_m"].values.copy()
    orig_uz = nodes["UZ_m"].values.copy()

    # Apply scale factor to displacements [m]
    nodes["UX_m"] = orig_ux * scale
    nodes["UY_m"] = orig_uy * scale
    nodes["UZ_m"] = orig_uz * scale

    grid = build_pipe_surface(data, deformed=True)

    # Restore original displacements
    nodes["UX_m"] = orig_ux
    nodes["UY_m"] = orig_uy
    nodes["UZ_m"] = orig_uz

    return grid


# ── Plotter Builder ───────────────────────────────────────────────────────


def build_plotter(
    data: dict,
    fields: dict[str, np.ndarray],
    result_name: str,
    colormap: str,
    deform_scale: float,
    show_labels: bool,
    undeformed_grid: pv.StructuredGrid,
) -> pv.Plotter:
    """Build a complete off-screen PyVista plotter with current settings.

    Args:
        data: Data dictionary.
        fields: Scalar fields dictionary.
        result_name: Display name of active result (key into ``RESULTS``).
        colormap: Matplotlib colormap name.
        deform_scale: Deformation magnification factor [-].
        show_labels: Whether to add min/max value labels on the surface.
        undeformed_grid: Pre-built undeformed surface (static, added as
            wireframe overlay).

    Returns:
        Configured off-screen PyVista Plotter.
    """
    result_info = RESULTS[result_name]
    scalar_key = result_info["key"]
    unit = result_info["unit"]

    # Build deformed surface with scaled displacements
    deformed = _build_scaled_surface(data, deform_scale)

    # Assign all scalar fields to deformed mesh
    for key, values in fields.items():
        deformed.point_data[key] = values

    pl = pv.Plotter(off_screen=True)

    # Undeformed wireframe (grey, 50% opacity)
    pl.add_mesh(
        undeformed_grid,
        style="wireframe",
        color="grey",
        opacity=0.5,
    )

    # Deformed surface coloured by active scalar
    pl.add_mesh(
        deformed,
        scalars=scalar_key,
        cmap=colormap,
        scalar_bar_args={
            "title": f"{result_name} [{unit}]",
            "shadow": True,
        },
    )

    # Fault zone (translucent red box)
    meta = data["metadata"]
    if "fault_zone_x_min" in meta and "fault_zone_x_max" in meta:
        x_min = float(meta["fault_zone_x_min"])
        x_max = float(meta["fault_zone_x_max"])
        r = data["outer_R"]  # [m]
        box = pv.Box(bounds=[x_min, x_max, -3 * r, 3 * r, -3 * r, 3 * r])
        pl.add_mesh(box, color="red", opacity=0.15)

    # Yielded spring markers (red spheres below deformed pipe)
    vt_mask, _ = _compute_yield_masks(data)
    if np.any(vt_mask):
        nodes_df = data["nodes_df"]
        xs = (
            nodes_df["x_m"].values[vt_mask]
            + nodes_df["UX_m"].values[vt_mask] * deform_scale
        )
        ys = nodes_df["UY_m"].values[vt_mask] * deform_scale - data["outer_R"]
        zs = nodes_df["UZ_m"].values[vt_mask] * deform_scale
        pl.add_mesh(
            pv.PolyData(np.column_stack([xs, ys, zs])),
            color="red",
            point_size=10,
            render_points_as_spheres=True,
        )

    # Min/max value labels on the pipe surface
    if show_labels:
        _add_minmax_labels(pl, data, fields, scalar_key, unit)

    pl.view_isometric()
    return pl


def _add_minmax_labels(
    pl: pv.Plotter,
    data: dict,
    fields: dict[str, np.ndarray],
    scalar_key: str,
    unit: str,
) -> None:
    """Add min/max value point labels to the plotter.

    Args:
        pl: PyVista plotter to add labels to.
        data: Data dictionary.
        fields: Scalar fields dictionary.
        scalar_key: Internal key of the active scalar.
        unit: Unit string for label formatting.
    """
    vals = fields[scalar_key]
    nc = N_CIRC + 1
    x_nodes = data["nodes_df"]["x_m"].values
    n_nodes = len(x_nodes)
    r = data["outer_R"]  # [m]

    # Per-node extreme value (max absolute across circumference)
    node_vals = np.array(
        [
            vals[i * nc : (i + 1) * nc][np.argmax(np.abs(vals[i * nc : (i + 1) * nc]))]
            for i in range(n_nodes)
        ]
    )

    max_idx = int(np.argmax(np.abs(node_vals)))
    min_idx = int(np.argmin(np.abs(node_vals)))

    max_val = node_vals[max_idx]
    min_val = node_vals[min_idx]
    max_x = x_nodes[max_idx]  # [m]
    min_x = x_nodes[min_idx]  # [m]

    pl.add_point_labels(
        np.array([[max_x, r * 1.5, 0]]),
        [f"Max: {max_val:.3e} {unit} @ x={max_x:.2f}m"],
        font_size=12,
        text_color="white",
        shape_color="black",
        shape="rounded_rect",
        shape_opacity=0.8,
    )
    if max_idx != min_idx:
        pl.add_point_labels(
            np.array([[min_x, -r * 1.5, 0]]),
            [f"Min: {min_val:.3e} {unit} @ x={min_x:.2f}m"],
            font_size=12,
            text_color="white",
            shape_color="darkblue",
            shape="rounded_rect",
            shape_opacity=0.8,
        )


# ── Summary Table ─────────────────────────────────────────────────────────


def make_summary_df(
    data: dict,
    fields: dict[str, np.ndarray],
    result_name: str,
) -> pd.DataFrame:
    """Create a one-row summary DataFrame for the active result.

    Args:
        data: Data dictionary.
        fields: Scalar fields dictionary.
        result_name: Display name of active result.

    Returns:
        DataFrame with columns: Result, Max Value, Max Location (x),
        Min Value, Min Location (x), Units.
    """
    info = RESULTS[result_name]
    scalar_key = info["key"]
    unit = info["unit"]
    vals = fields[scalar_key]
    nc = N_CIRC + 1
    x_nodes = data["nodes_df"]["x_m"].values
    n_nodes = len(x_nodes)

    # Per-node extreme values (max absolute across circumference)
    node_vals = np.array(
        [
            vals[i * nc : (i + 1) * nc][np.argmax(np.abs(vals[i * nc : (i + 1) * nc]))]
            for i in range(n_nodes)
        ]
    )

    max_idx = int(np.argmax(node_vals))
    min_idx = int(np.argmin(node_vals))

    return pd.DataFrame(
        {
            "Result": [result_name],
            "Max Value": [f"{node_vals[max_idx]:.4e}"],
            "Max Location (x)": [f"{x_nodes[max_idx]:.2f} m"],
            "Min Value": [f"{node_vals[min_idx]:.4e}"],
            "Min Location (x)": [f"{x_nodes[min_idx]:.2f} m"],
            "Units": [unit],
        }
    )


# ── Camera Presets ────────────────────────────────────────────────────────


def _camera_presets(mesh: pv.DataSet) -> dict[str, dict]:
    """Compute standard engineering camera views from mesh bounds.

    Each view is a dict with ``position``, ``focalPoint``, and ``viewUp``
    keys suitable for assignment to ``pn.pane.VTK.camera``.

    Args:
        mesh: PyVista mesh to compute views for.

    Returns:
        Dictionary mapping view names to camera parameter dicts.
    """
    bounds = mesh.bounds  # [xmin, xmax, ymin, ymax, zmin, zmax]
    cx = (bounds[0] + bounds[1]) / 2.0
    cy = (bounds[2] + bounds[3]) / 2.0
    cz = (bounds[4] + bounds[5]) / 2.0
    center = [cx, cy, cz]

    # Diagonal length of bounding box [m]
    diag = np.sqrt(
        (bounds[1] - bounds[0]) ** 2
        + (bounds[3] - bounds[2]) ** 2
        + (bounds[5] - bounds[4]) ** 2
    )
    d = diag * 1.5  # camera distance [m]
    s = 0.577  # 1/sqrt(3) for isometric

    return {
        "Top": {
            "position": [cx, cy, cz + d],
            "focalPoint": center,
            "viewUp": [0, 1, 0],
        },
        "Bot": {
            "position": [cx, cy, cz - d],
            "focalPoint": center,
            "viewUp": [0, 1, 0],
        },
        "Right": {
            "position": [cx, cy + d, cz],
            "focalPoint": center,
            "viewUp": [0, 0, 1],
        },
        "Left": {
            "position": [cx, cy - d, cz],
            "focalPoint": center,
            "viewUp": [0, 0, 1],
        },
        "Front": {
            "position": [cx + d, cy, cz],
            "focalPoint": center,
            "viewUp": [0, 0, 1],
        },
        "Back": {
            "position": [cx - d, cy, cz],
            "focalPoint": center,
            "viewUp": [0, 0, 1],
        },
        "Iso": {
            "position": [cx + d * s, cy + d * s, cz + d * s],
            "focalPoint": center,
            "viewUp": [0, 0, 1],
        },
    }


# ── Panel App Builder ────────────────────────────────────────────────────


def create_app(data: dict) -> pn.Column:
    """Build the Panel application layout with all widgets and callbacks.

    Args:
        data: Data dictionary with loaded and derived quantities.

    Returns:
        Panel Column layout ready to be served via ``pn.serve()``.
    """
    fields = compute_scalar_fields(data)
    undeformed_grid = build_pipe_surface(data, deformed=False)

    # ── Widgets ───────────────────────────────────────────────────────

    result_select = pn.widgets.Select(
        name="Result",
        options=list(RESULTS.keys()),
        value="Longitudinal Stress",
        width=220,
    )

    colormap_select = pn.widgets.Select(
        name="Colormap",
        options=["jet", "coolwarm", "viridis", "plasma", "RdBu"],
        value="jet",
        width=120,
    )

    deform_scale = pn.widgets.FloatSlider(
        name="Deform Scale",
        start=1.0,
        end=100.0,
        step=1.0,
        value=10.0,
        width=200,
    )

    btn_top = pn.widgets.Button(name="Top", button_type="default", width=60)
    btn_bot = pn.widgets.Button(name="Bot", button_type="default", width=60)
    btn_left = pn.widgets.Button(name="Left", button_type="default", width=60)
    btn_right = pn.widgets.Button(name="Right", button_type="default", width=60)
    btn_front = pn.widgets.Button(name="Front", button_type="default", width=60)
    btn_back = pn.widgets.Button(name="Back", button_type="default", width=60)
    btn_iso = pn.widgets.Button(name="Iso", button_type="primary", width=60)

    btn_screenshot = pn.widgets.Button(
        name="Save Screenshot", button_type="success", width=140
    )
    btn_labels = pn.widgets.Toggle(name="Show Min/Max Labels", value=False, width=160)

    btn_export_html = pn.widgets.Button(
        name="Export HTML", button_type="warning", width=120
    )
    status_text = pn.pane.Markdown(
        "_Ready_", width=300, styles={"font-size": "12px", "color": "#666"}
    )

    # ── Initial scene ─────────────────────────────────────────────────

    initial_pl = build_plotter(
        data,
        fields,
        result_select.value,
        colormap_select.value,
        deform_scale.value,
        btn_labels.value,
        undeformed_grid,
    )

    cameras = _camera_presets(undeformed_grid)

    vtk_pane = pn.pane.VTK(
        initial_pl.ren_win,
        sizing_mode="stretch_both",
        height=600,
    )

    summary_df = make_summary_df(data, fields, result_select.value)
    summary_table = pn.pane.DataFrame(summary_df, sizing_mode="stretch_width")

    # ── Rebuild callback ──────────────────────────────────────────────

    def _rebuild_scene(*_events: object) -> None:
        """Rebuild the 3D scene with current widget values."""
        result_name = result_select.value
        scalar_key = RESULTS[result_name]["key"]

        # Guard: skip all-zero or all-NaN fields
        vals = fields[scalar_key]
        if np.all(vals == 0) or np.all(np.isnan(vals)):
            print(f"  WARNING: {scalar_key} has no data — skipping")
            return

        pl = build_plotter(
            data,
            fields,
            result_name,
            colormap_select.value,
            deform_scale.value,
            btn_labels.value,
            undeformed_grid,
        )
        vtk_pane.object = pl.ren_win
        summary_table.object = make_summary_df(data, fields, result_name)

    result_select.param.watch(_rebuild_scene, "value")
    colormap_select.param.watch(_rebuild_scene, "value")
    deform_scale.param.watch(_rebuild_scene, "value_throttled")
    btn_labels.param.watch(_rebuild_scene, "value")

    # ── Camera button callbacks ───────────────────────────────────────

    def _make_camera_cb(view_name: str):
        """Return a click callback that sets the camera to *view_name*."""

        def _cb(_event: object) -> None:
            vtk_pane.camera = cameras[view_name]

        return _cb

    btn_top.on_click(_make_camera_cb("Top"))
    btn_bot.on_click(_make_camera_cb("Bot"))
    btn_left.on_click(_make_camera_cb("Left"))
    btn_right.on_click(_make_camera_cb("Right"))
    btn_front.on_click(_make_camera_cb("Front"))
    btn_back.on_click(_make_camera_cb("Back"))
    btn_iso.on_click(_make_camera_cb("Iso"))

    # ── Screenshot callback ───────────────────────────────────────────

    def _save_screenshot(_event: object) -> None:
        """Save a PNG screenshot to the results folder."""
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        rk = RESULTS[result_select.value]["key"]
        fname = f"viewer_screenshot_{rk}_{ts}.png"
        fpath = Path(data["results_dir"]) / fname
        pl = build_plotter(
            data,
            fields,
            result_select.value,
            colormap_select.value,
            deform_scale.value,
            btn_labels.value,
            undeformed_grid,
        )
        pl.screenshot(str(fpath))
        print(f"  Screenshot saved: {fpath}")

    btn_screenshot.on_click(_save_screenshot)

    # ── HTML export callback ──────────────────────────────────────────

    def _export_html(_event: object) -> None:
        """Export the 3D scene as a standalone HTML file."""
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        rname = result_select.value.replace(" ", "_").lower()
        fname = f"pipeline_viewer_{rname}_{ts}.html"
        fpath = str(Path(data["results_dir"]) / fname)

        print()
        print("=" * 60)
        print("Exporting standalone HTML...")
        print(f"Result: {result_select.value}")
        print(f"File: {fpath}")

        try:
            pl = build_plotter(
                data,
                fields,
                result_select.value,
                colormap_select.value,
                deform_scale.value,
                btn_labels.value,
                undeformed_grid,
            )
            pl.export_html(fpath)
            size_mb = os.path.getsize(fpath) / 1e6
            print(f"Size: {size_mb:.1f} MB")
            print("Export complete. Open in any browser — no Python required.")
            print("=" * 60)
            status_text.object = f"Exported: {os.path.basename(fpath)}"
            if size_mb > 50:
                print("WARNING: File is large — may be slow to open in browser")
        except Exception as exc:
            print(f"ERROR: PyVista HTML export failed: {exc}")
            print("Try installing: pip install nest_asyncio")
            print("=" * 60)
            status_text.object = "Export failed — see terminal"

    btn_export_html.on_click(_export_html)

    # ── Assembly ──────────────────────────────────────────────────────

    toolbar = pn.Row(
        result_select,
        colormap_select,
        deform_scale,
        pn.Spacer(width=20),
        btn_top,
        btn_bot,
        btn_left,
        btn_right,
        btn_front,
        btn_back,
        btn_iso,
        pn.Spacer(width=20),
        btn_screenshot,
        btn_labels,
        btn_export_html,
        status_text,
        sizing_mode="stretch_width",
        margin=(5, 10),
    )

    # Store a reference so the Panel fallback export can access the layout
    app_layout = pn.Column(
        pn.pane.Markdown("## OpenSees Pipeline Viewer", margin=(5, 10)),
        toolbar,
        vtk_pane,
        summary_table,
        sizing_mode="stretch_both",
    )

    return app_layout


# ── Entry Point ───────────────────────────────────────────────────────────


def main() -> None:
    """Entry point for the web-based pipeline viewer."""
    if len(sys.argv) > 1:
        results_dir = sys.argv[1].strip().strip("'\"")
    else:
        results_dir = input("Drag results folder here, then press Enter: ").strip()
        results_dir = results_dir.strip("'\"")

    if not os.path.isdir(results_dir):
        print(f"ERROR: Directory not found: {results_dir}")
        sys.exit(1)

    data = load_results(results_dir)
    data = compute_derived(data)

    # Print the same terminal summary as the desktop viewer
    fields = compute_scalar_fields(data)
    print_summary(data, fields)

    print()
    print("=" * 60)
    print("OpenSees Pipeline Viewer — Web UI")
    print("=" * 60)
    print(f"Results loaded from: {data['results_dir']}")
    print("Opening browser at: http://localhost:5007")
    print("Press Ctrl+C to stop the server.")
    print("=" * 60)

    pn.extension("vtk")

    app = create_app(data)
    pn.serve(app, port=5007, show=True, title="OpenSees Pipeline Viewer")


if __name__ == "__main__":
    main()

# Required packages: pip install panel pyvista trame trame-vtk trame-vuetify
