cmake_minimum_required(VERSION 3.16)
project(CFDSolver LANGUAGES CXX)

set(PYBIND11_FINDPYTHON ON)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
add_subdirectory(external/pybind11)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

enable_testing()

add_library(cfd_linalg
    src/linalg/gaussian_elimination.cpp
    src/linalg/matrix.cpp
    src/linalg/vector.cpp
)

add_library(cfd_core
    src/core/velocity_field.cpp
    src/core/pressure_field.cpp
    src/core/boundary_conditions.cpp
    src/core/grid.cpp
    src/core/simulator.cpp
)

pybind11_add_module(cfdsolver_py
    src/python/bindings.cpp
)

target_link_libraries(cfdsolver_py PRIVATE cfd_core)
target_include_directories(cfdsolver_py PRIVATE include)


target_include_directories(cfd_linalg PUBLIC include)
target_include_directories(cfd_core PUBLIC include)

# Unit tests
add_executable(test_gaussian_elimination
    tests/unit_tests/test_gaussian_elimination.cpp
)
add_executable(test_velocity_field 
    tests/unit_tests/test_velocity_field.cpp
)
add_executable(test_pressure_field 
    tests/unit_tests/test_pressure_field.cpp
)
add_executable(test_boundary_conditions 
    tests/unit_tests/test_boundary_conditions.cpp
)
add_executable(test_grid 
    tests/unit_tests/test_grid.cpp
)

# Integration tests
add_executable(test_invariants
    tests/integration_tests/test_invariants.cpp
)
add_executable(test_dynamics
    tests/integration_tests/test_dynamics.cpp
)

target_link_libraries(cfd_core PUBLIC cfd_linalg)
target_link_libraries(test_gaussian_elimination PRIVATE cfd_linalg)
target_link_libraries(test_velocity_field PRIVATE cfd_core)
target_link_libraries(test_pressure_field PRIVATE cfd_core)
target_link_libraries(test_boundary_conditions PRIVATE cfd_core)
target_link_libraries(test_grid PRIVATE cfd_core)
target_link_libraries(test_invariants PRIVATE cfd_core)
target_link_libraries(test_dynamics PRIVATE cfd_core)


# Unit tests
add_test(
    NAME test_gaussian_elimination
    COMMAND test_gaussian_elimination
)
add_test(
    NAME test_velocity_field
    COMMAND test_velocity_field
)
add_test(
    NAME test_pressure_field
    COMMAND test_pressure_field
)
add_test(
    NAME test_boundary_conditions
    COMMAND test_boundary_conditions
)
add_test(
    NAME test_grid
    COMMAND test_grid
)

# Integration tests
add_test(
    NAME test_invariants
    COMMAND test_invariants
)

add_test(
    NAME test_dynamics
    COMMAND test_dynamics
)

# Benchmarks

add_executable(lid_driven_cavity_benchmark
    benchmarks/lid_driven_cavity.cpp
)

target_link_libraries(lid_driven_cavity_benchmark PRIVATE cfd_core)