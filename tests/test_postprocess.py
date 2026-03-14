"""Tests for post-processing module.

These tests are placeholders for when the postprocess module is fully
implemented. They document the expected interface and behavior.
"""

from __future__ import annotations

import pytest

from opensees_pipeline.postprocess import (
    compute_bending_stress,
    extract_bending_moments,
    extract_displacements,
)


class TestExtractDisplacements:
    """Tests for displacement extraction (stub)."""

    def test_not_implemented(self) -> None:
        """Function raises NotImplementedError until implemented."""
        with pytest.raises(NotImplementedError):
            extract_displacements()


class TestExtractBendingMoments:
    """Tests for moment extraction (stub)."""

    def test_not_implemented(self) -> None:
        """Function raises NotImplementedError until implemented."""
        with pytest.raises(NotImplementedError):
            extract_bending_moments()


class TestComputeBendingStress:
    """Tests for stress computation (stub)."""

    def test_not_implemented(self) -> None:
        """Function raises NotImplementedError until implemented."""
        with pytest.raises(NotImplementedError):
            compute_bending_stress()
