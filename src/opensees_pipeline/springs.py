"""ASCE ALA 2001 soil spring backbone curve calculations.

This module implements the bilinear elastic-perfectly-plastic force-displacement
backbone curves for buried pipe-soil interaction springs following ASCE ALA (2001)
"Guidelines for the Design of Buried Steel Pipe", Sections 7.2–7.5.

Four spring directions are provided:
    - Axial (t-z): resistance to pipe movement along its axis
    - Lateral horizontal (p-y): resistance to pipe movement perpendicular to axis
    - Vertical bearing (q-z): resistance to downward pipe movement
    - Uplift (q-z upward): resistance to upward pipe movement

All functions return *per-unit-length* values. These must be multiplied by the
tributary length at each node before being assigned to OpenSees zero-length
elements. This scaling is typically done in the model assembly step.

References:
    ASCE ALA (2001). Guidelines for the Design of Buried Steel Pipe.
    American Lifelines Alliance, ASCE.
"""

from __future__ import annotations

import math
from typing import Literal


def ala2001_axial_spring(
    pipe_od: float,
    burial_depth: float,
    soil_unit_weight: float,
    friction_angle: float,
    interface_friction_factor: float,
    adhesion: float = 0.0,
    unit_system: Literal["SI", "Imperial"] = "SI",
) -> tuple[float, float]:
    """Compute axial (t-z) spring parameters per ALA 2001 Section 7.2.

    The axial soil spring resists relative movement between the pipe and soil
    along the pipe axis. The ultimate axial force per unit length is derived
    from the frictional resistance on the pipe-soil interface.

    ALA 2001 Eq. (referenced in Section 7.2):
        tu = pi * D * (1/2 * gamma * H * (1 + K0) * tan(delta) + c_a)

    where:
        - D = pipe outside diameter
        - gamma = soil unit weight
        - H = depth to pipe centerline
        - K0 = coefficient of lateral earth pressure at rest = 1 - sin(phi)
        - delta = interface friction angle = f * phi
        - c_a = soil-pipe adhesion for cohesive soils
        - f = interface friction factor (typically 0.5–1.0)

    The yield displacement zu is taken as a small value (typically 3–5 mm for
    dense sand, 8–10 mm for loose sand in SI units).

    Args:
        pipe_od: Pipe outside diameter [m | ft].
        burial_depth: Depth from ground surface to pipe centerline [m | ft].
        soil_unit_weight: Total unit weight of soil [kN/m^3 | pcf].
        friction_angle: Soil friction angle [degrees].
        interface_friction_factor: Ratio of interface friction angle to soil
            friction angle, f = delta/phi. Typically 0.5–1.0 depending on
            pipe coating. ALA 2001 Table 7-1.
        adhesion: Soil-pipe adhesion for cohesive soils [kPa | psf].
            Default 0.0 for granular soils.
        unit_system: "SI" for (m, kN, kPa) or "Imperial" for (ft, kips, ksf).

    Returns:
        Tuple of (tu, zu) where:
            tu: Ultimate axial force per unit length [kN/m | kip/ft].
            zu: Yield displacement [m | ft].

    Note:
        Results are PER UNIT LENGTH. Multiply by tributary node length before
        assigning to OpenSees zero-length elements.
    """
    phi_rad = math.radians(friction_angle)  # friction angle [rad]
    delta_rad = interface_friction_factor * phi_rad  # interface friction angle [rad]
    k0 = 1.0 - math.sin(phi_rad)  # at-rest earth pressure coefficient [-]

    # ALA 2001 Section 7.2: tu = pi * D * (0.5 * gamma * H * (1 + K0) * tan(delta) + c_a)
    tu = (
        math.pi
        * pipe_od  # D [m | ft]
        * (
            0.5
            * soil_unit_weight  # gamma [kN/m^3 | pcf]
            * burial_depth  # H [m | ft]
            * (1.0 + k0)  # (1 + K0) [-]
            * math.tan(delta_rad)  # tan(delta) [-]
            + adhesion  # c_a [kPa | psf]
        )
    )

    # Yield displacement: ALA 2001 recommends 3-5 mm for dense sand, 8-10 mm for loose
    if unit_system == "SI":
        zu = 0.003  # 3 mm [m] — conservative for dense sand
    else:
        zu = 0.01  # ~3 mm [ft] (0.12 in)

    return tu, zu


def ala2001_lateral_spring(
    pipe_od: float,
    burial_depth: float,
    soil_unit_weight: float,
    friction_angle: float,
    cohesion: float = 0.0,
    unit_system: Literal["SI", "Imperial"] = "SI",
) -> tuple[float, float]:
    """Compute lateral horizontal (p-y) spring parameters per ALA 2001 Section 7.3.

    The lateral soil spring resists pipe movement perpendicular to its axis in
    the horizontal plane. The ultimate lateral force per unit length depends on
    the passive earth pressure acting on the pipe cross-section.

    ALA 2001 Section 7.3:
        pu = N_ch * c * D + N_qh * gamma * H * D

    where:
        - N_ch = horizontal bearing capacity factor for cohesion (function of H/D)
        - N_qh = horizontal bearing capacity factor for overburden (function of H/D, phi)
        - c = soil cohesion
        - D = pipe outside diameter
        - gamma = soil unit weight
        - H = depth to pipe centerline

    The yield displacement is typically 0.04*(H + D/2) for loose to medium sand
    and 0.02*(H + D/2) for dense sand.

    Args:
        pipe_od: Pipe outside diameter [m | ft].
        burial_depth: Depth from ground surface to pipe centerline [m | ft].
        soil_unit_weight: Total unit weight of soil [kN/m^3 | pcf].
        friction_angle: Soil friction angle [degrees].
        cohesion: Soil cohesion [kPa | psf]. Default 0.0 for granular soils.
        unit_system: "SI" for (m, kN, kPa) or "Imperial" for (ft, kips, ksf).

    Returns:
        Tuple of (pu, yu) where:
            pu: Ultimate lateral force per unit length [kN/m | kip/ft].
            yu: Yield displacement [m | ft].

    Note:
        Results are PER UNIT LENGTH. Multiply by tributary node length before
        assigning to OpenSees zero-length elements.
    """
    h_over_d = burial_depth / pipe_od  # embedment ratio [-]

    # N_ch: horizontal bearing capacity factor for cohesion
    # ALA 2001 Figure 7.3 — varies with H/D, capped at 9
    n_ch = min(6.752 * h_over_d**0.345, 9.0)  # [-] empirical fit

    # N_qh: horizontal bearing capacity factor for overburden
    # ALA 2001 Figure 7.4 — depends on H/D and phi
    # Simplified closed-form approximation per Hansen (1961) as referenced in ALA
    phi_rad = math.radians(friction_angle)  # [rad]
    # Approximate N_qh using ALA recommended values
    # For granular soils, N_qh increases with H/D and phi
    kp = (1.0 + math.sin(phi_rad)) / (1.0 - math.sin(phi_rad))  # passive earth pressure [-]
    n_qh = h_over_d * kp  # simplified approximation [-]
    # Cap N_qh per ALA 2001 recommendations (varies with phi)
    n_qh_max = {20: 2.4, 25: 3.9, 30: 6.7, 35: 10.8, 40: 17.5, 45: 28.0}
    # Interpolate cap from nearest phi values
    phi_deg = friction_angle
    cap_keys = sorted(n_qh_max.keys())
    if phi_deg <= cap_keys[0]:
        cap = n_qh_max[cap_keys[0]]
    elif phi_deg >= cap_keys[-1]:
        cap = n_qh_max[cap_keys[-1]]
    else:
        for i in range(len(cap_keys) - 1):
            if cap_keys[i] <= phi_deg <= cap_keys[i + 1]:
                frac = (phi_deg - cap_keys[i]) / (cap_keys[i + 1] - cap_keys[i])
                cap = n_qh_max[cap_keys[i]] + frac * (
                    n_qh_max[cap_keys[i + 1]] - n_qh_max[cap_keys[i]]
                )
                break
    n_qh = min(n_qh, cap)

    # ALA 2001 Eq: pu = N_ch * c * D + N_qh * gamma * H * D
    pu = (
        n_ch * cohesion * pipe_od  # cohesive component [kN/m | kip/ft]
        + n_qh * soil_unit_weight * burial_depth * pipe_od  # frictional component
    )

    # Yield displacement: 0.04 * (H + D/2) for medium-density sand
    yu = 0.04 * (burial_depth + pipe_od / 2.0)  # [m | ft]

    return pu, yu


def ala2001_bearing_spring(
    pipe_od: float,
    burial_depth: float,
    soil_unit_weight: float,
    friction_angle: float,
    cohesion: float = 0.0,
    unit_system: Literal["SI", "Imperial"] = "SI",
) -> tuple[float, float]:
    """Compute vertical bearing (q-z downward) spring parameters per ALA 2001 Section 7.4.

    The vertical bearing spring resists downward pipe movement into the soil.
    The ultimate bearing capacity is based on the general bearing capacity
    equation adapted for a buried pipe.

    ALA 2001 Section 7.4:
        qu = N_c * c * D + N_q * gamma * H * D + 0.5 * N_gamma * gamma * D^2

    where:
        - N_c, N_q, N_gamma = bearing capacity factors (functions of phi)
        - c = soil cohesion
        - D = pipe outside diameter
        - gamma = soil unit weight
        - H = depth to pipe centerline

    The yield displacement is typically 0.10 * D for granular soils and
    0.20 * D for cohesive soils.

    Args:
        pipe_od: Pipe outside diameter [m | ft].
        burial_depth: Depth from ground surface to pipe centerline [m | ft].
        soil_unit_weight: Total unit weight of soil [kN/m^3 | pcf].
        friction_angle: Soil friction angle [degrees].
        cohesion: Soil cohesion [kPa | psf]. Default 0.0 for granular soils.
        unit_system: "SI" for (m, kN, kPa) or "Imperial" for (ft, kips, ksf).

    Returns:
        Tuple of (qu, zq) where:
            qu: Ultimate bearing force per unit length [kN/m | kip/ft].
            zq: Yield displacement [m | ft].

    Note:
        Results are PER UNIT LENGTH. Multiply by tributary node length before
        assigning to OpenSees zero-length elements.
    """
    phi_rad = math.radians(friction_angle)  # [rad]

    # Bearing capacity factors (Terzaghi/Meyerhof as referenced in ALA 2001)
    # N_q = e^(pi*tan(phi)) * tan^2(45 + phi/2)
    n_q = math.exp(math.pi * math.tan(phi_rad)) * math.tan(
        math.radians(45.0 + friction_angle / 2.0)
    ) ** 2  # [-]
    # N_c = (N_q - 1) / tan(phi) — L'Hôpital for phi=0 gives N_c = 5.14
    if friction_angle > 0.1:
        n_c = (n_q - 1.0) / math.tan(phi_rad)  # [-]
    else:
        n_c = 5.14  # [-] undrained case
    # N_gamma = 2 * (N_q + 1) * tan(phi) — Vesic approximation
    n_gamma = 2.0 * (n_q + 1.0) * math.tan(phi_rad)  # [-]

    # ALA 2001 Eq: qu = N_c * c * D + N_q * gamma * H * D + 0.5 * N_gamma * gamma * D^2
    qu = (
        n_c * cohesion * pipe_od  # cohesion term [kN/m | kip/ft]
        + n_q * soil_unit_weight * burial_depth * pipe_od  # overburden term
        + 0.5 * n_gamma * soil_unit_weight * pipe_od**2  # self-weight term
    )

    # Yield displacement: 0.10 * D for granular soils
    zq = 0.10 * pipe_od  # [m | ft]

    return qu, zq


def ala2001_uplift_spring(
    pipe_od: float,
    burial_depth: float,
    soil_unit_weight: float,
    friction_angle: float,
    cohesion: float = 0.0,
    unit_system: Literal["SI", "Imperial"] = "SI",
) -> tuple[float, float]:
    """Compute vertical uplift (q-z upward) spring parameters per ALA 2001 Section 7.5.

    The uplift spring resists upward pipe movement out of the soil. The ultimate
    uplift resistance is governed by the weight of the soil wedge above the pipe
    plus any cohesive resistance along the failure surface.

    ALA 2001 Section 7.5:
        Qu = N_cv * c * D + N_qv * gamma * H * D

    where:
        - N_cv = vertical uplift factor for cohesion (function of H/D)
        - N_qv = vertical uplift factor for overburden (function of H/D, phi)
        - c = soil cohesion
        - D = pipe outside diameter
        - gamma = soil unit weight
        - H = depth to pipe centerline

    The total uplift resistance must not exceed the weight of the soil column
    above the pipe: Qu_max = gamma * H * D (for granular soils).

    The yield displacement is typically 0.01*H to 0.02*H for dense sand
    and 0.02*H to 0.05*H for loose sand.

    Args:
        pipe_od: Pipe outside diameter [m | ft].
        burial_depth: Depth from ground surface to pipe centerline [m | ft].
        soil_unit_weight: Total unit weight of soil [kN/m^3 | pcf].
        friction_angle: Soil friction angle [degrees].
        cohesion: Soil cohesion [kPa | psf]. Default 0.0 for granular soils.
        unit_system: "SI" for (m, kN, kPa) or "Imperial" for (ft, kips, ksf).

    Returns:
        Tuple of (qu_up, zqu) where:
            qu_up: Ultimate uplift force per unit length [kN/m | kip/ft].
            zqu: Yield displacement [m | ft].

    Note:
        Results are PER UNIT LENGTH. Multiply by tributary node length before
        assigning to OpenSees zero-length elements.
    """
    h_over_d = burial_depth / pipe_od  # embedment ratio [-]
    phi_rad = math.radians(friction_angle)  # [rad]

    # N_cv: vertical uplift factor for cohesion
    # ALA 2001 Section 7.5 — similar to N_ch but for vertical direction
    n_cv = min(2.0 * h_over_d, 10.0)  # [-] capped at 10

    # N_qv: vertical uplift factor for overburden
    # ALA 2001 Section 7.5 — function of H/D and phi
    # Simplified: N_qv = phi/44 * (H/D) for H/D <= 11, capped
    n_qv = (friction_angle / 44.0) * h_over_d  # [-]
    # Cap N_qv per ALA 2001 recommendations
    n_qv_max = {20: 1.6, 25: 3.2, 30: 5.5, 35: 8.5, 40: 11.5, 45: 15.0}
    phi_deg = friction_angle
    cap_keys = sorted(n_qv_max.keys())
    if phi_deg <= cap_keys[0]:
        cap = n_qv_max[cap_keys[0]]
    elif phi_deg >= cap_keys[-1]:
        cap = n_qv_max[cap_keys[-1]]
    else:
        for i in range(len(cap_keys) - 1):
            if cap_keys[i] <= phi_deg <= cap_keys[i + 1]:
                frac = (phi_deg - cap_keys[i]) / (cap_keys[i + 1] - cap_keys[i])
                cap = n_qv_max[cap_keys[i]] + frac * (
                    n_qv_max[cap_keys[i + 1]] - n_qv_max[cap_keys[i]]
                )
                break
    n_qv = min(n_qv, cap)

    # ALA 2001 Eq: Qu = N_cv * c * D + N_qv * gamma * H * D
    qu_up = (
        n_cv * cohesion * pipe_od  # cohesive component [kN/m | kip/ft]
        + n_qv * soil_unit_weight * burial_depth * pipe_od  # frictional component
    )

    # Cap at weight of soil column above pipe: gamma * H * D
    qu_max = soil_unit_weight * burial_depth * pipe_od  # [kN/m | kip/ft]
    qu_up = min(qu_up, qu_max)

    # Yield displacement: 0.02 * H for medium-density sand
    zqu = 0.02 * burial_depth  # [m | ft]

    return qu_up, zqu
