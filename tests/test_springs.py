"""Tests for ALA 2001 soil spring calculations."""

from __future__ import annotations

import math

import pytest

from opensees_pipeline.springs import (
    ala2001_axial_spring,
    ala2001_bearing_spring,
    ala2001_lateral_spring,
    ala2001_uplift_spring,
)


# Common test parameters — typical buried steel pipe in medium-dense sand
PIPE_OD = 0.3048  # 12-inch pipe [m]
BURIAL_DEPTH = 1.5  # depth to centerline [m]
SOIL_GAMMA = 18.0  # unit weight [kN/m^3]
PHI = 35.0  # friction angle [degrees]


class TestAxialSpring:
    """Tests for ala2001_axial_spring."""

    def test_positive_force(self) -> None:
        """Ultimate axial force must be positive."""
        tu, zu = ala2001_axial_spring(
            PIPE_OD, BURIAL_DEPTH, SOIL_GAMMA, PHI, 0.7
        )
        assert tu > 0.0

    def test_positive_displacement(self) -> None:
        """Yield displacement must be positive."""
        tu, zu = ala2001_axial_spring(
            PIPE_OD, BURIAL_DEPTH, SOIL_GAMMA, PHI, 0.7
        )
        assert zu > 0.0

    def test_force_increases_with_depth(self) -> None:
        """Deeper burial produces higher axial resistance."""
        tu_shallow, _ = ala2001_axial_spring(PIPE_OD, 1.0, SOIL_GAMMA, PHI, 0.7)
        tu_deep, _ = ala2001_axial_spring(PIPE_OD, 3.0, SOIL_GAMMA, PHI, 0.7)
        assert tu_deep > tu_shallow

    def test_force_increases_with_friction(self) -> None:
        """Higher interface friction factor produces higher axial resistance."""
        tu_low, _ = ala2001_axial_spring(PIPE_OD, BURIAL_DEPTH, SOIL_GAMMA, PHI, 0.5)
        tu_high, _ = ala2001_axial_spring(PIPE_OD, BURIAL_DEPTH, SOIL_GAMMA, PHI, 1.0)
        assert tu_high > tu_low

    def test_zero_friction_angle(self) -> None:
        """Zero friction angle with adhesion still produces resistance."""
        tu, zu = ala2001_axial_spring(
            PIPE_OD, BURIAL_DEPTH, SOIL_GAMMA, 0.0, 0.7, adhesion=10.0
        )
        assert tu > 0.0

    def test_imperial_units(self) -> None:
        """Imperial unit system returns different yield displacement."""
        _, zu_si = ala2001_axial_spring(
            PIPE_OD, BURIAL_DEPTH, SOIL_GAMMA, PHI, 0.7, unit_system="SI"
        )
        _, zu_imp = ala2001_axial_spring(
            1.0, 5.0, 120.0, PHI, 0.7, unit_system="Imperial"
        )
        assert zu_si != zu_imp


class TestLateralSpring:
    """Tests for ala2001_lateral_spring."""

    def test_positive_force(self) -> None:
        """Ultimate lateral force must be positive."""
        pu, yu = ala2001_lateral_spring(PIPE_OD, BURIAL_DEPTH, SOIL_GAMMA, PHI)
        assert pu > 0.0

    def test_positive_displacement(self) -> None:
        """Yield displacement must be positive."""
        pu, yu = ala2001_lateral_spring(PIPE_OD, BURIAL_DEPTH, SOIL_GAMMA, PHI)
        assert yu > 0.0

    def test_force_increases_with_depth(self) -> None:
        """Deeper burial produces higher lateral resistance."""
        pu_shallow, _ = ala2001_lateral_spring(PIPE_OD, 1.0, SOIL_GAMMA, PHI)
        pu_deep, _ = ala2001_lateral_spring(PIPE_OD, 3.0, SOIL_GAMMA, PHI)
        assert pu_deep > pu_shallow

    def test_cohesion_contribution(self) -> None:
        """Adding cohesion increases lateral resistance."""
        pu_no_c, _ = ala2001_lateral_spring(PIPE_OD, BURIAL_DEPTH, SOIL_GAMMA, PHI)
        pu_with_c, _ = ala2001_lateral_spring(
            PIPE_OD, BURIAL_DEPTH, SOIL_GAMMA, PHI, cohesion=20.0
        )
        assert pu_with_c > pu_no_c


class TestBearingSpring:
    """Tests for ala2001_bearing_spring."""

    def test_positive_force(self) -> None:
        """Ultimate bearing force must be positive."""
        qu, zq = ala2001_bearing_spring(PIPE_OD, BURIAL_DEPTH, SOIL_GAMMA, PHI)
        assert qu > 0.0

    def test_positive_displacement(self) -> None:
        """Yield displacement must be positive."""
        qu, zq = ala2001_bearing_spring(PIPE_OD, BURIAL_DEPTH, SOIL_GAMMA, PHI)
        assert zq > 0.0

    def test_bearing_exceeds_uplift(self) -> None:
        """Bearing capacity should exceed uplift capacity for same conditions."""
        qu, _ = ala2001_bearing_spring(PIPE_OD, BURIAL_DEPTH, SOIL_GAMMA, PHI)
        qu_up, _ = ala2001_uplift_spring(PIPE_OD, BURIAL_DEPTH, SOIL_GAMMA, PHI)
        assert qu > qu_up

    def test_yield_displacement_proportional_to_diameter(self) -> None:
        """Bearing yield displacement is 0.10 * D."""
        _, zq = ala2001_bearing_spring(PIPE_OD, BURIAL_DEPTH, SOIL_GAMMA, PHI)
        assert zq == pytest.approx(0.10 * PIPE_OD)


class TestUpliftSpring:
    """Tests for ala2001_uplift_spring."""

    def test_positive_force(self) -> None:
        """Ultimate uplift force must be positive."""
        qu_up, zqu = ala2001_uplift_spring(PIPE_OD, BURIAL_DEPTH, SOIL_GAMMA, PHI)
        assert qu_up > 0.0

    def test_positive_displacement(self) -> None:
        """Yield displacement must be positive."""
        qu_up, zqu = ala2001_uplift_spring(PIPE_OD, BURIAL_DEPTH, SOIL_GAMMA, PHI)
        assert zqu > 0.0

    def test_capped_at_soil_weight(self) -> None:
        """Uplift resistance must not exceed soil weight above pipe."""
        qu_up, _ = ala2001_uplift_spring(PIPE_OD, BURIAL_DEPTH, SOIL_GAMMA, PHI)
        qu_max = SOIL_GAMMA * BURIAL_DEPTH * PIPE_OD  # [kN/m]
        assert qu_up <= qu_max + 1e-10

    def test_force_increases_with_depth(self) -> None:
        """Deeper burial produces higher uplift resistance."""
        qu_shallow, _ = ala2001_uplift_spring(PIPE_OD, 0.5, SOIL_GAMMA, PHI)
        qu_deep, _ = ala2001_uplift_spring(PIPE_OD, 3.0, SOIL_GAMMA, PHI)
        assert qu_deep > qu_shallow
