from typing import Any


SETUP_DEFAULTS: dict[str, dict[str, Any]] = {
    "Lid-driven cavity": {
        "grid_w": 100,
        "grid_h": 100,
        "length": 1.0,
        "rho": 1.225,
        "reynolds": 100.0,
        "lid_velocity": 1.0,
    },

    "Pressure-driven channel": {
        "grid_w": 120,
        "grid_h": 60,
        "length": 1.0,
        "rho": 1.225,
        "viscosity": 0.1,
        "pressure_difference": 10.0,
    },

    "Airfoil flow": {
        "grid_w": 250,
        "grid_h": 200,
        "length": 1.0,
        "rho": 1.225,
        "viscosity": 0.1,
        "naca_code": "2412",
        "chord": 0.5,
        "freestream_velocity": 1.0,
        "angle_of_attack_deg": 5.0,
    },
}


PARAMS_BY_SETUP: dict[str, list[str]] = {
    "Lid-driven cavity": [
        "grid_w",
        "grid_h",
        "length",
        "rho",
        "reynolds",
        "lid_velocity",
    ],

    "Pressure-driven channel": [
        "grid_w",
        "grid_h",
        "length",
        "rho",
        "viscosity",
        "pressure_difference",
    ],

    "Airfoil flow": [
        "grid_w",
        "grid_h",
        "length",
        "rho",
        "viscosity",
        "naca_code",
        "chord",
        "freestream_velocity",
        "angle_of_attack_deg",
    ],
}


PARAM_LABELS: dict[str, str] = {
    "grid_w": "Grid width",
    "grid_h": "Grid height",
    "length": "Domain length",
    "rho": "Density",
    "reynolds": "Reynolds number",
    "lid_velocity": "Lid velocity",
    "viscosity": "Viscosity",
    "pressure_difference": "Pressure difference",
    "naca_code": "NACA code",
    "chord": "Chord",
    "freestream_velocity": "Freestream velocity",
    "angle_of_attack_deg": "Angle of attack",
}


INT_PARAMS: set[str] = {
    "grid_w",
    "grid_h",
}


TEXT_PARAMS: set[str] = {
    "naca_code",
}


# Per-parameter (min, max) overrides for spin boxes.
# Float params default to (1e-9, 1e9); int params default to (1, 10_000).
PARAM_RANGES: dict[str, tuple[float, float]] = {
    "angle_of_attack_deg": (-90.0, 90.0),
    "pressure_difference": (-1e6, 1e6),
}


def decimals_for_param(name: str) -> int:
    if name in {"length", "chord"}:
        return 3

    if name in {"rho", "viscosity"}:
        return 4

    if name in {"lid_velocity", "freestream_velocity"}:
        return 2

    if name in {"pressure_difference", "angle_of_attack_deg"}:
        return 2

    if name == "reynolds":
        return 1

    return 3


def step_for_param(name: str) -> float:
    if name == "reynolds":
        return 10.0

    if name == "pressure_difference":
        return 1.0

    if name == "angle_of_attack_deg":
        return 1.0

    if name in {"length", "chord"}:
        return 0.01

    if name in {"rho", "viscosity"}:
        return 0.001

    if name in {"lid_velocity", "freestream_velocity"}:
        return 0.1

    return 0.1