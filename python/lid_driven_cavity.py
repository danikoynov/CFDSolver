import os
import sys
from pathlib import Path

root = Path(__file__).resolve().parent.parent
build_dir = root / "build"

if hasattr(os, "add_dll_directory"):
    os.add_dll_directory(r"C:\msys64\ucrt64\bin")

sys.path.insert(0, str(build_dir))

import cfdsolver_py
from visualizer import Visualizer

def impose_boundary_conditions(sim: cfdsolver_py.Simulator, lid_velocity: float):
    grid = sim.grid()
    bc = grid.boundary_conditions()

    w = grid.width()
    h = grid.height()
    
    for i in range(w):
        bc.set_cell_type(i, 0, cfdsolver_py.CellType.SOLID)
        bc.set_cell_type(i, h - 1, cfdsolver_py.CellType.SOLID)

        bc.prescribe_u_value(i + 1, 0, 0)
        bc.prescribe_u_value(i + 1, h - 1, 0)

        bc.prescribe_v_value(i, 1, 0)
        bc.prescribe_v_value(i, h - 1, 0)

    for i in range(2, w - 1):
        bc.prescribe_u_value(i, h - 2, lid_velocity)

    for j in range(h):
        bc.set_cell_type(0, j, cfdsolver_py.CellType.SOLID)
        bc.set_cell_type(w - 1, j, cfdsolver_py.CellType.SOLID)

        bc.prescribe_v_value(0, j + 1, 0)
        bc.prescribe_v_value(w - 1, j + 1, 0)

        bc.prescribe_u_value(1, j, 0)
        bc.prescribe_u_value(w - 1, j, 0)
        

if __name__ == "__main__":
    
    length = 1.0 # Side of cavity
    lid_velocity = 1.0 # Lid velocity

    Re = 100 # Re = length * lid_velocity / viscosity
    viscosity = length * lid_velocity / Re
    
    grid_w = 20
    grid_h = 20
    cell_size = 40  # reduced from 40 so the larger grid fits better

    dx = length / grid_w
    
    sim = cfdsolver_py.Simulator(
        grid_w,
        grid_h,
        dx,
        0.0,
        1.225,
        False,
        viscosity
    )
    
    impose_boundary_conditions(sim, lid_velocity)
    
    visualizer = Visualizer(sim, grid_w, grid_h, cell_size)
     
    visualizer.run()
    
    
