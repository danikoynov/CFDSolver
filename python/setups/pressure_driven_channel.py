import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt

here = Path(__file__).resolve()
python_dir = here.parent.parent          # .../python
root = here.parent.parent.parent         # project root
build_dir = root / "build"

if hasattr(os, "add_dll_directory"):
    os.add_dll_directory(r"C:\msys64\ucrt64\bin")

sys.path.insert(0, str(build_dir))
sys.path.insert(0, str(python_dir))

import cfdsolver_py
from visualizer import Visualizer


# --------------------------------------------------
# Choose mode here:
# "visualizer" -> run pygame visualizer
# "profile"    -> run simulation and plot velocity profile
# --------------------------------------------------
RUN_MODE = "profile"


def impose_boundary_conditions(sim: cfdsolver_py.Simulator, pd: float):
    grid = sim.grid()
    bc = grid.boundary_conditions()

    w = grid.width()
    h = grid.height()

    bc.prescribe_u_value(0, 0, 0)
    bc.prescribe_u_value(0, h - 1, 0)

    for i in range(w):
        bc.set_cell_type(i, 0, cfdsolver_py.CellType.SOLID)
        bc.set_cell_type(i, h - 1, cfdsolver_py.CellType.SOLID)

        bc.prescribe_u_value(i + 1, 0, 0)
        bc.prescribe_u_value(i + 1, h - 1, 0)

        bc.prescribe_v_value(i, 1, 0)
        bc.prescribe_v_value(i, h - 1, 0)

    for j in range(1, h - 1):
        bc.prescribe_p_value(-1, j, pd)
        bc.prescribe_p_value(w, j, 0)


def cell_center_velocity(vf, x, y):
    u = 0.5 * (vf.get_u(x, y) + vf.get_u(x + 1, y))
    v = 0.5 * (vf.get_v(x, y) + vf.get_v(x, y + 1))
    return u, v


import numpy as np
import matplotlib.pyplot as plt


def plot_velocity_profile(sim, grid_w, grid_h, nsteps=500):
    # advance simulation
    for _ in range(nsteps):
        sim.tick()

    grid = sim.grid()
    vf = grid.velocity()

    # vertical centerline: plot u(y) at x = mid
    x_mid = grid_w // 2

    y_vals = []
    u_vals = []

    for y in range(grid_h):
        u, _ = cell_center_velocity(vf, x_mid, y)

        # normalize y to [0, 1]
        y_center = (y + 0.5) / grid_h

        y_vals.append(y_center)
        u_vals.append(u)

    # fit a parabola u(y) = a y^2 + b y + c
    coeffs = np.polyfit(y_vals, u_vals, 2)
    a, b, c = coeffs

    # smooth y values for plotting fitted curve
    y_fit = np.linspace(min(y_vals), max(y_vals), 200)
    u_fit = a * y_fit**2 + b * y_fit + c

    plt.figure(figsize=(6, 6))
    plt.plot(u_vals, y_vals, "o", markersize=3, label="numerical")
    plt.plot(u_fit, y_fit, "-", linewidth=2, label="quadratic fit")

    plt.xlabel("u velocity")
    plt.ylabel("y")
    plt.title(f"Vertical centerline velocity profile after {nsteps} steps")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    print(f"Fitted parabola: u(y) = {a:.6f} y^2 + {b:.6f} y + {c:.6f}")

if __name__ == "__main__":

    length = 1.0
    viscosity = 0.1

    grid_w = 50
    grid_h = 50
    cell_size = 12

    pd = 10.0
    dx = length / grid_w

    sim = cfdsolver_py.Simulator(
        grid_w,
        grid_h,
        dx,
        1.225,
        False,
        viscosity
    )

    impose_boundary_conditions(sim, pd)

    if RUN_MODE == "visualizer":
        visualizer = Visualizer(sim, grid_w, grid_h, cell_size)
        visualizer.run()

    elif RUN_MODE == "profile":
        plot_velocity_profile(sim, grid_w, grid_h, nsteps=500)

    else:
        raise ValueError(
            f"Unknown RUN_MODE = {RUN_MODE!r}. "
            f'Use "visualizer" or "profile".'
        )