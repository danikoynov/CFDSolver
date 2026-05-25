# FlowGrid — Improvement Roadmap

## Performance — Critical

- [ ] **Warm-start conjugate gradient** (`src/linalg/conjugate_gradient.cpp`)
  Start each tick's CG solve from the previous tick's pressure solution instead of x=0.
  For low-Re flows where pressure changes slowly between ticks, this can cut iteration
  counts by 50–80%.

- [ ] **Cache the Poisson matrix between ticks** (`src/core/simulator.cpp` → `project()`)
  The matrix is rebuilt from scratch every tick even though the stencil only changes
  when the grid topology changes (i.e., on setup, not during simulation). Store it as a
  member and rebuild only when `configure()` is called.

## Performance — High

- [ ] **Replace element-wise CG convergence check with L2 norm** (`src/linalg/conjugate_gradient.cpp` → `within_tolerance()`)
  The current check iterates over every element to find the max absolute value. A single
  dot-product norm (`r·r < tol² * b·b`) is cheaper and is the standard stopping criterion.

- [ ] **Bulk numpy export in C++ bindings** (`src/python/bindings.cpp`)
  `pressure_field()` and `velocity_fields()` in `SimulationController` call `p.get(i, j)`
  and `vf.get_u/v` per cell in Python — O(N) pybind11 round-trips per frame. Add
  `pressure_array()` and `velocity_arrays()` methods to the C++ simulator that return
  `py::array_t<float>` directly, reducing render overhead to a single call.

## Performance — Medium

- [ ] **Reuse VelocityField buffers in `advect()` and `apply_viscosity()`** (`src/core/simulator.cpp`)
  Both functions allocate a new `VelocityField` every tick. Store a scratch buffer as a
  member and swap in-place to eliminate per-tick heap allocation.

- [ ] **Multithreaded conjugate gradient** (`src/linalg/conjugate_gradient.cpp`)
  The dot products and SpMV in the CG loop are trivially parallelisable with OpenMP or
  `std::execution::par_unseq`. Noted in the README as future work.

## Correctness

- [ ] **`AirfoilFlowSetup` should inherit `BaseSetup` in pybind11 bindings** (`src/python/bindings.cpp`)
  All other setup classes use `py::class_<..., cfd::setups::BaseSetup>` as their base.
  `AirfoilFlowSetup` is registered without this base, breaking Python-side polymorphism
  and making it impossible to accept any `BaseSetup&` parameter from Python.

- [ ] **Validate viscosity solver stability** (`src/core/simulator.cpp` → `apply_viscosity()`)
  The explicit diffusion step is conditionally stable only when `dt ≤ h²/(2ν)`. The
  timestep code already computes a diffusive limit, but it should be checked that this
  limit is always active and not overridden by the CFL bound alone at high viscosity.

## Features

- [ ] **Field data export from GUI** (`python/gui/`)
  Add an "Export" button that dumps the current pressure/velocity arrays to `.npy` or
  CSV. Useful for post-processing and regression testing.

- [ ] **Validation benchmarks** (`tests/`)
  Add Ghia et al. (1982) reference data for Re=400 and Re=1000 lid-driven cavity and
  assert that the steady-state u-velocity profile along the vertical centreline matches
  within a tolerance. Currently the only automated test is a basic smoke test.

- [ ] **Inflow / outflow boundary conditions**
  The channel setup uses a pressure difference but has no true inflow BC. A Dirichlet
  velocity inlet and a zero-gradient Neumann outlet would enable more realistic channel
  and external flow setups.

- [ ] **Pressure-driven channel: compute viscosity from Reynolds number**
  Like the lid-driven cavity, the channel setup should let the user specify Re and derive
  viscosity automatically rather than requiring the user to enter it manually.

## UX / GUI

- [ ] **Configurable render frame rate** (`python/gui/main_window.py`)
  The timer is hardcoded at 33 ms (~30 fps). An fps slider in the SIMULATION section
  (e.g. 10–60 fps) would let users trade visual smoothness for raw throughput.

- [ ] **Elapsed simulation time display** (`python/gui/main_window.py`)
  Show wall-clock time alongside the iteration counter so the user can gauge solver speed.

- [ ] **Setup change resets iteration counter visually**
  After `apply_configuration()` the counter resets to 0 but the reset isn't immediately
  obvious. A brief colour flash on the label would make the transition clear.
