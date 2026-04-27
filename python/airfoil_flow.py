import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt

here = Path(__file__).resolve()
python_dir = here.parent          # .../python
root = here.parent.parent         # project root
build_dir = root / "build"

if hasattr(os, "add_dll_directory"):
    os.add_dll_directory(r"C:\msys64\ucrt64\bin")

sys.path.insert(0, str(build_dir))
sys.path.insert(0, str(python_dir))

import cfdsolver_py
from visualizer import Visualizer
from airfoil_generator.naca_airfoil import Airfoil

def impose_boundary_conditions(
    sim: cfdsolver_py.Simulator,
    airfoil: Airfoil,
    v_air: float,
    dx: float,
):
    grid = sim.grid()
    bc = grid.boundary_conditions()
    vf = grid.velocity()

    w = grid.width()
    h = grid.height()

    domain_width = w * dx
    domain_height = h * dx

    # --------------------------------------------------
    # Initialize the whole flow field to freestream
    # This avoids the "still fluid + sudden inlet" pressure jump.
    # --------------------------------------------------
    for j in range(h):
        for i in range(w + 1):
            vf.set_u(i, j, v_air)

    for j in range(h + 1):
        for i in range(w):
            vf.set_v(i, j, 0.0)

    # --------------------------------------------------
    # Airfoil placement: centered, AoA = 0
    # --------------------------------------------------
    x_le = 0.5 * domain_width - 0.5 * airfoil.chord
    y_chord = 0.5 * domain_height

    def cell_center(i: int, j: int):
        x = (i + 0.5) * dx
        y = (j + 0.5) * dx
        return x, y

    def is_inside_airfoil(i: int, j: int) -> bool:
        x, y = cell_center(i, j)

        X = (x - x_le) / airfoil.chord

        if X < 0.0 or X > 1.0:
            return False

        y_c = airfoil.mean_camber_line(X)
        y_t = airfoil.half_thickness(X)

        y_airfoil_center = y_chord + y_c

        return (y_airfoil_center - y_t) <= y <= (y_airfoil_center + y_t)

    def make_solid_cell(i: int, j: int):
        bc.set_cell_type(i, j, cfdsolver_py.CellType.SOLID)

        # Set actual velocity field to zero on the solid cell faces
        vf.set_u(i, j, 0.0)
        vf.set_u(i + 1, j, 0.0)

        vf.set_v(i, j, 0.0)
        vf.set_v(i, j + 1, 0.0)

        # Also prescribe those values as boundary conditions
        bc.prescribe_u_value(i, j, 0.0)
        bc.prescribe_u_value(i + 1, j, 0.0)

        bc.prescribe_v_value(i, j, 0.0)
        bc.prescribe_v_value(i, j + 1, 0.0)

    # --------------------------------------------------
    # Left inlet velocity
    # --------------------------------------------------
    for j in range(h):
        bc.prescribe_u_value(0, j, v_air)

    # --------------------------------------------------
    # Top and bottom far-field
    # --------------------------------------------------
    for i in range(w):
        # no normal flow through top/bottom
        bc.prescribe_v_value(i, 0, 0.0)
        bc.prescribe_v_value(i, h, 0.0)

        # horizontal freestream along top/bottom
        bc.prescribe_u_value(i, 0, v_air)
        bc.prescribe_u_value(i, h - 1, v_air)

    # --------------------------------------------------
    # Right pressure outlet
    # Uses your ghost-cell pressure binding.
    # --------------------------------------------------
    for j in range(h):
        bc.prescribe_p_value(w, j, 0.0)
        bc.prescribe_p_value(-1, j, 0.0)

    # --------------------------------------------------
    # Airfoil solid cells
    # --------------------------------------------------
    solid_count = 0

    for j in range(h):
        for i in range(w):
            if is_inside_airfoil(i, j):
                make_solid_cell(i, j)
                solid_count += 1

    print(f"Marked {solid_count} airfoil cells as SOLID")
    
def cell_center_velocity(vf, x, y):
    u = 0.5 * (vf.get_u(x, y) + vf.get_u(x + 1, y))
    v = 0.5 * (vf.get_v(x, y) + vf.get_v(x, y + 1))
    return u, v


if __name__ == "__main__":

    length = 1.0
    viscosity = 0.1

    grid_w = 110
    grid_h = 40
    cell_size = 10

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

    airfoil = Airfoil("2412", 0.5)
    v_air = 1 # m/s
    
    impose_boundary_conditions(sim, airfoil, v_air, dx)


    scale_arrows = True
    visualizer = Visualizer(sim, grid_w, grid_h, cell_size, scale_arrows)
    visualizer.run()