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
# "profile"    -> run simulation and plot centerline profiles
# --------------------------------------------------
RUN_MODE = "visualizer"

def cell_center_velocity(vf, x, y):
    u = 0.5 * (vf.get_u(x, y) + vf.get_u(x + 1, y))
    v = 0.5 * (vf.get_v(x, y) + vf.get_v(x, y + 1))
    return u, v


def plot_velocity_profiles(sim, grid_w, grid_h, nsteps=1000):
    for _ in range(nsteps):
        sim.tick()

    grid = sim.grid()
    vf = grid.velocity()

    x_mid = grid_w // 2
    y_mid = grid_h // 2

    # vertical centerline: u(y) at x = mid
    y_vals = []
    u_vals = []
    for y in range(grid_h):
        u, _ = cell_center_velocity(vf, x_mid, y)
        y_center = (y + 0.5) / grid_h
        y_vals.append(y_center)
        u_vals.append(u)

    # horizontal centerline: v(x) at y = mid
    x_vals = []
    v_vals = []
    for x in range(grid_w):
        _, v = cell_center_velocity(vf, x, y_mid)
        x_center = (x + 0.5) / grid_w
        x_vals.append(x_center)
        v_vals.append(v)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

    ax1.plot(u_vals, y_vals, marker="o", markersize=3)
    ax1.set_xlabel("u velocity")
    ax1.set_ylabel("y")
    ax1.set_title(f"Vertical centerline u(y) after {nsteps} steps")
    ax1.grid(True)

    ax2.plot(x_vals, v_vals, marker="o", markersize=3)
    ax2.set_xlabel("x")
    ax2.set_ylabel("v velocity")
    ax2.set_title(f"Horizontal centerline v(x) after {nsteps} steps")
    ax2.grid(True)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    length = 1.0
    lid_velocity = 1.0

    Re = 100
    viscosity = length * lid_velocity / Re

    grid_w = 100
    grid_h = 100
    cell_size = 8

    dx = length / grid_w

    sim = cfdsolver_py.Simulator(
        grid_w,
        grid_h,
        dx,
        1.225,
        False,
        viscosity
    )
    setup = cfdsolver_py.LidDrivenCavity(lid_velocity)
    setup.impose_boundary_conditions(sim.grid())
    
    if RUN_MODE == "visualizer":
        scale_arrows = False
        visualizer = Visualizer(sim, grid_w, grid_h, cell_size, scale_arrows)
        visualizer.run()

    elif RUN_MODE == "profile":
        plot_velocity_profiles(sim, grid_w, grid_h, nsteps=1000)

    else:
        raise ValueError(
            f"Unknown RUN_MODE = {RUN_MODE!r}. "
            f'Use "visualizer" or "profile".'
        )