"""Tests for mesh generation module."""

from __future__ import annotations

import numpy as np
import pytest

from opensees_pipeline.mesh import (
    biased_spacing,
    compute_tributary_lengths,
    generate_mesh_from_kps,
)


class TestBiasedSpacing:
    """Tests for biased_spacing function."""

    def test_uniform_spacing(self) -> None:
        """Uniform spacing (bias_ratio=1.0) produces evenly spaced nodes."""
        coords = biased_spacing(4, 10.0, bias_ratio=1.0)
        expected = np.array([0.0, 2.5, 5.0, 7.5, 10.0])
        np.testing.assert_allclose(coords, expected)

    def test_node_count(self) -> None:
        """Number of nodes is n_elem + 1."""
        coords = biased_spacing(6, 12.0, bias_ratio=2.0)
        assert len(coords) == 7

    def test_endpoints(self) -> None:
        """First node at 0, last node at length."""
        coords = biased_spacing(5, 15.0, bias_ratio=3.0)
        assert coords[0] == pytest.approx(0.0)
        assert coords[-1] == pytest.approx(15.0)

    def test_bias_fine_at_start(self) -> None:
        """With fine_end='start', first element is smallest."""
        coords = biased_spacing(5, 10.0, bias_ratio=3.0, fine_end="start")
        element_lengths = np.diff(coords)
        assert element_lengths[0] < element_lengths[-1]

    def test_bias_fine_at_end(self) -> None:
        """With fine_end='end', last element is smallest."""
        coords = biased_spacing(5, 10.0, bias_ratio=3.0, fine_end="end")
        element_lengths = np.diff(coords)
        assert element_lengths[-1] < element_lengths[0]

    def test_single_element(self) -> None:
        """Single element produces two nodes at endpoints."""
        coords = biased_spacing(1, 5.0)
        np.testing.assert_allclose(coords, [0.0, 5.0])

    def test_invalid_n_elem(self) -> None:
        """Raises ValueError for n_elem < 1."""
        with pytest.raises(ValueError, match="n_elem"):
            biased_spacing(0, 10.0)

    def test_invalid_length(self) -> None:
        """Raises ValueError for non-positive length."""
        with pytest.raises(ValueError, match="length"):
            biased_spacing(5, -1.0)

    def test_invalid_bias_ratio(self) -> None:
        """Raises ValueError for bias_ratio < 1."""
        with pytest.raises(ValueError, match="bias_ratio"):
            biased_spacing(5, 10.0, bias_ratio=0.5)


class TestComputeTributaryLengths:
    """Tests for compute_tributary_lengths function."""

    def test_uniform_mesh(self) -> None:
        """Uniform mesh: interior nodes have elem_length, end nodes have half."""
        x = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
        trib = compute_tributary_lengths(x)
        expected = np.array([0.5, 1.0, 1.0, 1.0, 0.5])
        np.testing.assert_allclose(trib, expected)

    def test_sum_equals_total_length(self) -> None:
        """Sum of tributary lengths equals total pipe length."""
        x = np.array([0.0, 1.0, 2.5, 5.0, 8.0])
        trib = compute_tributary_lengths(x)
        assert trib.sum() == pytest.approx(x[-1] - x[0])

    def test_two_nodes(self) -> None:
        """Two nodes: each gets half the element length."""
        x = np.array([0.0, 6.0])
        trib = compute_tributary_lengths(x)
        np.testing.assert_allclose(trib, [3.0, 3.0])

    def test_nonuniform_mesh(self) -> None:
        """Non-uniform mesh produces correct tributary lengths."""
        x = np.array([0.0, 1.0, 3.0, 6.0])
        trib = compute_tributary_lengths(x)
        # Node 0: 1.0/2 = 0.5
        # Node 1: (1.0 + 2.0)/2 = 1.5
        # Node 2: (2.0 + 3.0)/2 = 2.5
        # Node 3: 3.0/2 = 1.5
        expected = np.array([0.5, 1.5, 2.5, 1.5])
        np.testing.assert_allclose(trib, expected)

    def test_too_few_nodes(self) -> None:
        """Raises ValueError for fewer than 2 nodes."""
        with pytest.raises(ValueError, match="At least 2 nodes"):
            compute_tributary_lengths(np.array([0.0]))


class TestGenerateMeshFromKps:
    """Tests for generate_mesh_from_kps function."""

    def test_two_kps_straight(self) -> None:
        """Two key points produce a straight mesh."""
        kps = [(0.0, 0.0, 0.0), (10.0, 0.0, 0.0)]
        nodes, elements = generate_mesh_from_kps(kps, default_elem_length=2.5)
        assert nodes.shape[0] == 5  # 4 elements + 1
        assert len(elements) == 4

    def test_node_endpoints(self) -> None:
        """First and last nodes match key points."""
        kps = [(0.0, 0.0, 0.0), (10.0, 5.0, 0.0)]
        nodes, _ = generate_mesh_from_kps(kps, default_elem_length=2.0)
        np.testing.assert_allclose(nodes[0], [0.0, 0.0, 0.0])
        np.testing.assert_allclose(nodes[-1], [10.0, 5.0, 0.0], atol=1e-10)

    def test_three_kps_no_duplicates(self) -> None:
        """Three key points: shared node is not duplicated."""
        kps = [(0.0, 0.0, 0.0), (5.0, 0.0, 0.0), (10.0, 0.0, 0.0)]
        nodes, elements = generate_mesh_from_kps(kps, default_elem_length=5.0)
        # Each segment: 1 element, so 2+1=3 nodes but shared = 3 nodes total
        # Actually 5/5=1 elem per seg, so 2 nodes per seg, minus 1 shared = 3
        assert nodes.shape[0] == 3
        assert len(elements) == 2

    def test_too_few_kps(self) -> None:
        """Raises ValueError for fewer than 2 key points."""
        with pytest.raises(ValueError, match="At least 2 key points"):
            generate_mesh_from_kps([(0.0, 0.0, 0.0)], default_elem_length=1.0)
