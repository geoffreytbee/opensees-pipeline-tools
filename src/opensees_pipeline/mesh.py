"""Mesh generation for buried pipeline models.

This module provides functions for generating finite element meshes along 3D
pipeline alignments. It handles node placement, element connectivity, biased
mesh refinement near features of interest, and tributary length computation
for soil spring scaling.

Tributary lengths are critical for buried pipeline analysis because soil spring
properties (from ALA 2001 or PRCI) are defined per unit length of pipe. When
assigning springs to discrete nodes in a finite element model, the per-unit-length
values must be scaled by the tributary length — the length of pipe "served" by
each node. For interior nodes, this is the average of the two adjacent element
lengths. For end nodes, it is half the adjacent element length.
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd


def biased_spacing(
    n_elem: int,
    length: float,
    bias_ratio: float = 1.0,
    fine_end: Literal["start", "end"] = "start",
) -> np.ndarray:
    """Generate geometrically biased element spacing along a segment.

    Creates a geometric progression of element lengths, with smaller elements
    concentrated at one end. This is useful for mesh refinement near features
    such as fault crossings, pipe bends, or boundary conditions where stress
    gradients are steep.

    The bias ratio is defined as the ratio of the largest element to the
    smallest element. A bias_ratio of 1.0 produces uniform spacing.

    Args:
        n_elem: Number of elements in the segment.
        length: Total length of the segment [m | ft].
        bias_ratio: Ratio of largest to smallest element length [-].
            Must be >= 1.0. Default 1.0 (uniform spacing).
        fine_end: Which end has the finest (smallest) elements.
            "start" puts small elements at x=0, "end" puts them at x=length.

    Returns:
        Array of node coordinates along the segment, shape (n_elem + 1,),
        starting at 0.0 and ending at length.

    Raises:
        ValueError: If n_elem < 1, length <= 0, or bias_ratio < 1.0.

    Example:
        >>> coords = biased_spacing(5, 10.0, bias_ratio=3.0, fine_end="start")
        >>> len(coords)
        6
        >>> coords[0], coords[-1]
        (0.0, 10.0)
    """
    if n_elem < 1:
        raise ValueError(f"n_elem must be >= 1, got {n_elem}")
    if length <= 0:
        raise ValueError(f"length must be > 0, got {length}")
    if bias_ratio < 1.0:
        raise ValueError(f"bias_ratio must be >= 1.0, got {bias_ratio}")

    if n_elem == 1 or abs(bias_ratio - 1.0) < 1e-10:
        # Uniform spacing
        return np.linspace(0.0, length, n_elem + 1)

    # Geometric progression: element i has length a * r^i
    # where r = bias_ratio^(1/(n-1)) and a = length * (1-r)/(1-r^n)
    r = bias_ratio ** (1.0 / (n_elem - 1))  # common ratio [-]
    a = length * (1.0 - r) / (1.0 - r**n_elem)  # first element length [m | ft]

    # Cumulative positions
    element_lengths = a * r ** np.arange(n_elem)  # [m | ft] per element
    coords = np.zeros(n_elem + 1)
    coords[1:] = np.cumsum(element_lengths)
    coords[-1] = length  # enforce exact endpoint

    if fine_end == "end":
        # Reverse so fine elements are at the far end
        coords = length - coords[::-1]

    return coords


def generate_mesh_from_kps(
    kp_coords: list[tuple[float, float, float]],
    default_elem_length: float,
    feature_zones: list[dict] | None = None,
) -> tuple[np.ndarray, list[tuple[int, int]]]:
    """Generate node coordinates and element connectivity from key points.

    Creates a 3D pipeline mesh by placing nodes along straight segments between
    key points, with optional mesh refinement in specified feature zones.

    Args:
        kp_coords: List of (x, y, z) key point coordinates defining the
            pipeline alignment [m | ft]. Must have at least 2 points.
        default_elem_length: Default target element length [m | ft].
            Actual lengths may vary slightly to fit segment lengths exactly.
        feature_zones: Optional list of feature zone dictionaries, each with:
            - "start_kp": int — index of starting key point
            - "end_kp": int — index of ending key point
            - "elem_length": float — target element length in zone [m | ft]
            - "bias_ratio": float — bias ratio for mesh grading (default 1.0)
            If None, uniform meshing with default_elem_length is used.

    Returns:
        Tuple of (nodes, elements) where:
            nodes: Array of shape (n_nodes, 3) with x, y, z coordinates.
            elements: List of (node_i, node_j) connectivity tuples (0-indexed).

    Raises:
        ValueError: If fewer than 2 key points are provided.
    """
    if len(kp_coords) < 2:
        raise ValueError("At least 2 key points required")

    kp_array = np.array(kp_coords, dtype=float)
    all_nodes: list[np.ndarray] = []
    all_elements: list[tuple[int, int]] = []
    node_offset = 0

    for seg_idx in range(len(kp_array) - 1):
        p1 = kp_array[seg_idx]
        p2 = kp_array[seg_idx + 1]
        seg_length = float(np.linalg.norm(p2 - p1))  # [m | ft]

        # Determine element length for this segment
        elem_len = default_elem_length  # [m | ft]
        bias = 1.0  # [-]
        if feature_zones:
            for zone in feature_zones:
                if zone["start_kp"] <= seg_idx < zone["end_kp"]:
                    elem_len = zone.get("elem_length", default_elem_length)
                    bias = zone.get("bias_ratio", 1.0)
                    break

        n_elem = max(1, round(seg_length / elem_len))

        # Generate parametric coordinates along segment
        if abs(bias - 1.0) < 1e-10:
            t_values = np.linspace(0.0, 1.0, n_elem + 1)
        else:
            local_coords = biased_spacing(n_elem, seg_length, bias)
            t_values = local_coords / seg_length

        # Interpolate 3D coordinates
        segment_nodes = p1[np.newaxis, :] + t_values[:, np.newaxis] * (p2 - p1)

        # Skip first node of subsequent segments to avoid duplicates
        if seg_idx > 0:
            segment_nodes = segment_nodes[1:]

        # Build connectivity
        for i in range(len(segment_nodes) - 1):
            ni = node_offset + i
            nj = node_offset + i + 1
            all_elements.append((ni, nj))

        all_nodes.append(segment_nodes)
        node_offset += len(segment_nodes)

    nodes = np.vstack(all_nodes)
    return nodes, all_elements


def read_alignment_csv(
    filepath: str | Path,
    unit_system: Literal["SI", "Imperial"] = "SI",
) -> list[tuple[float, float, float]]:
    """Read pipeline alignment coordinates from CSV or Excel file.

    Expects a file with columns named 'x', 'y', 'z' (case-insensitive)
    containing the pipeline alignment coordinates. Supports CSV (.csv)
    and Excel (.xlsx, .xls) formats.

    Args:
        filepath: Path to the alignment data file.
        unit_system: "SI" for meters or "Imperial" for feet. Used for
            labeling only — no conversion is performed. The coordinates
            in the file must already be in the target unit system.

    Returns:
        List of (x, y, z) coordinate tuples [m | ft].

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If required columns are missing.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Alignment file not found: {filepath}")

    # Read CSV or Excel based on extension
    if filepath.suffix.lower() in (".xlsx", ".xls"):
        df = pd.read_excel(filepath)
    else:
        df = pd.read_csv(filepath)

    # Normalize column names to lowercase
    df.columns = [c.strip().lower() for c in df.columns]

    required_cols = {"x", "y", "z"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}. Found: {list(df.columns)}")

    coords = list(zip(df["x"].values, df["y"].values, df["z"].values))
    return coords


def compute_tributary_lengths(x_nodes: np.ndarray) -> np.ndarray:
    """Compute tributary length at each node along a 1D mesh.

    Tributary length is the length of pipe "served" by each node. It determines
    how per-unit-length soil spring properties are scaled to per-node values
    for assignment to OpenSees zero-length elements.

    For interior nodes, the tributary length is the average of the two adjacent
    element lengths:
        L_trib[i] = (x[i] - x[i-1] + x[i+1] - x[i]) / 2

    For end nodes, the tributary length is half the adjacent element length:
        L_trib[0] = (x[1] - x[0]) / 2
        L_trib[-1] = (x[-1] - x[-2]) / 2

    This ensures that the sum of all tributary lengths equals the total pipe
    length, which is a necessary condition for the discrete spring model to
    be equivalent to the continuous spring model.

    Args:
        x_nodes: Array of node positions along the pipe axis [m | ft].
            Must be sorted in ascending order with at least 2 nodes.

    Returns:
        Array of tributary lengths at each node [m | ft], same length as x_nodes.

    Raises:
        ValueError: If fewer than 2 nodes are provided.

    Example:
        >>> import numpy as np
        >>> x = np.array([0.0, 1.0, 2.5, 5.0])
        >>> trib = compute_tributary_lengths(x)
        >>> trib  # [0.5, 1.25, 1.75, 1.25] — sums to 5.0 (total length)
    """
    if len(x_nodes) < 2:
        raise ValueError("At least 2 nodes required for tributary length calculation")

    n = len(x_nodes)
    trib = np.zeros(n)

    # Element lengths
    elem_lengths = np.diff(x_nodes)  # [m | ft]

    # End nodes: half of adjacent element
    trib[0] = elem_lengths[0] / 2.0  # [m | ft]
    trib[-1] = elem_lengths[-1] / 2.0  # [m | ft]

    # Interior nodes: average of adjacent elements
    for i in range(1, n - 1):
        trib[i] = (elem_lengths[i - 1] + elem_lengths[i]) / 2.0  # [m | ft]

    return trib
