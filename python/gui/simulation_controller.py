import os
import sys
from pathlib import Path
from typing import Any

import numpy as np

from gui.setup_schema import SETUP_DEFAULTS


# ----------------------------------------------------------------------
# Make C++ extension importable
# ----------------------------------------------------------------------

here = Path(__file__).resolve()
python_dir = here.parent.parent
root = python_dir.parent
build_dir = root / "build"

if hasattr(os, "add_dll_directory"):
    os.add_dll_directory(r"C:\msys64\ucrt64\bin")

sys.path.insert(0, str(build_dir))
sys.path.insert(0, str(python_dir))

import cfdsolver_py


def cell_center_velocity(vf, x: int, y: int) -> tuple[float, float]:
    u = 0.5 * (vf.get_u(x, y) + vf.get_u(x + 1, y))
    v = 0.5 * (vf.get_v(x, y) + vf.get_v(x, y + 1))
    return u, v


class SimulationController:
    def __init__(self):
        self.setup_name = "Lid-driven cavity"
        self.params: dict[str, Any] = dict(SETUP_DEFAULTS[self.setup_name])

        self.grid_w = 0
        self.grid_h = 0
        self.length = 1.0
        self.dx = 1.0
        self.rho = 1.225

        self.sim = None
        self.setup = None

        self.step_count = 0
        self._solid_cache: np.ndarray | None = None

        self.create_simulation()

    # ------------------------------------------------------------------
    # Setup/configuration
    # ------------------------------------------------------------------

    def available_setups(self) -> list[str]:
        return list(SETUP_DEFAULTS.keys())

    def default_params_for(self, setup_name: str) -> dict[str, Any]:
        if setup_name not in SETUP_DEFAULTS:
            raise ValueError(f"Unknown setup: {setup_name}")

        return dict(SETUP_DEFAULTS[setup_name])

    def configure(self, setup_name: str, params: dict[str, Any]) -> None:
        if setup_name not in SETUP_DEFAULTS:
            raise ValueError(f"Unknown setup: {setup_name}")

        self.setup_name = setup_name
        self.params = dict(params)

        self.create_simulation()

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def create_simulation(self) -> None:
        self.grid_w = int(self.params["grid_w"])
        self.grid_h = int(self.params["grid_h"])
        self.length = float(self.params["length"])

        if self.grid_w <= 0 or self.grid_h <= 0:
            raise ValueError("Grid width and height must be positive.")

        if self.length <= 0.0:
            raise ValueError("Domain length must be positive.")

        self.dx = self.length / self.grid_w
        self.rho = float(self.params["rho"])

        if self.rho <= 0.0:
            raise ValueError("Density must be positive.")

        viscosity = self._compute_viscosity()

        self.sim = cfdsolver_py.Simulator(
            self.grid_w,
            self.grid_h,
            self.dx,
            self.rho,
            False,
            viscosity,
        )

        self._apply_selected_setup()

        self.step_count = 0
        self._solid_cache = None

        self._jj, self._ii = np.mgrid[0:self.grid_h, 0:self.grid_w]
        self._j_flat = self._jj.ravel().tolist()
        self._i_flat = self._ii.ravel().tolist()

    def reset(self) -> None:
        self.create_simulation()

    def step(self, steps: int = 1) -> None:
        for _ in range(steps):
            self.sim.tick()
            self.step_count += 1

    # ------------------------------------------------------------------
    # Setup internals
    # ------------------------------------------------------------------

    def _compute_viscosity(self) -> float:
        if self.setup_name == "Lid-driven cavity":
            lid_velocity = float(self.params["lid_velocity"])
            reynolds = float(self.params["reynolds"])

            if reynolds <= 0.0:
                raise ValueError("Reynolds number must be positive.")

            return self.length * lid_velocity / reynolds

        if self.setup_name == "Pressure-driven channel":
            viscosity = float(self.params["viscosity"])

            if viscosity <= 0.0:
                raise ValueError("Viscosity must be positive.")

            return viscosity

        if self.setup_name == "Airfoil flow":
            viscosity = float(self.params["viscosity"])

            if viscosity <= 0.0:
                raise ValueError("Viscosity must be positive.")

            return viscosity

        raise ValueError(f"Unknown setup: {self.setup_name}")

    def _apply_selected_setup(self) -> None:
        if self.setup_name == "Lid-driven cavity":
            self.setup = cfdsolver_py.LidDrivenCavity(
                float(self.params["lid_velocity"])
            )

            self.setup.impose_boundary_conditions(self.sim.grid())
            return

        if self.setup_name == "Pressure-driven channel":
            self.setup = cfdsolver_py.PressureDrivenChannel(
                float(self.params["pressure_difference"])
            )

            self.setup.impose_boundary_conditions(self.sim.grid())
            return

        if self.setup_name == "Airfoil flow":
            naca_code = str(self.params["naca_code"])
            chord = float(self.params["chord"])
            freestream_velocity = float(self.params["freestream_velocity"])
            angle_of_attack_deg = float(self.params["angle_of_attack_deg"])

            if len(naca_code) != 4 or not naca_code.isdigit():
                raise ValueError("NACA code must be exactly 4 digits.")

            if chord <= 0.0:
                raise ValueError("Chord must be positive.")

            if chord >= self.length:
                raise ValueError("Chord should be smaller than the domain length.")

            if freestream_velocity < 0.0:
                raise ValueError("Freestream velocity cannot be negative.")

            self.setup = cfdsolver_py.AirfoilFlowSetup(
                naca_code,
                chord,
                freestream_velocity,
                self.dx,
                angle_of_attack_deg,
            )

            self.setup.impose_boundary_conditions(self.sim.grid())
            return

        raise ValueError(f"Unknown setup: {self.setup_name}")

    # ------------------------------------------------------------------
    # Field accessors
    # ------------------------------------------------------------------

    def _flat_to_grid(self, values) -> np.ndarray:
        return np.array(values, dtype=np.float32).reshape(
            self.grid_h,
            self.grid_w,
        )

    def pressure_field(self) -> np.ndarray:
        p = self.sim.grid().pressure()

        vals = [
            p.get(i, j)
            for j, i in zip(self._j_flat, self._i_flat)
        ]

        return self._flat_to_grid(vals)

    def velocity_fields(self) -> tuple[np.ndarray, np.ndarray]:
        vf = self.sim.grid().velocity()

        pairs = [
            cell_center_velocity(vf, i, j)
            for j, i in zip(self._j_flat, self._i_flat)
        ]

        u_flat, v_flat = zip(*pairs)

        u = np.array(u_flat, dtype=np.float32).reshape(
            self.grid_h,
            self.grid_w,
        )

        v = np.array(v_flat, dtype=np.float32).reshape(
            self.grid_h,
            self.grid_w,
        )

        return u, v

    def speed_field(self) -> np.ndarray:
        u, v = self.velocity_fields()
        return np.hypot(u, v)

    def solid_mask(self) -> np.ndarray:
        if self._solid_cache is not None:
            return self._solid_cache

        bc = self.sim.grid().boundary_conditions()

        vals = [
            bc.type(i, j) == cfdsolver_py.CellType.SOLID
            for j, i in zip(self._j_flat, self._i_flat)
        ]

        self._solid_cache = np.array(vals, dtype=bool).reshape(
            self.grid_h,
            self.grid_w,
        )

        return self._solid_cache