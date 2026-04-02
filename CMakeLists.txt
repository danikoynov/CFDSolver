cmake_minimum_required(VERSION 3.16)
project(CFDSolver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

enable_testing()

add_library(cfd_linalg
    src/linalg/gaussian_elimination.cpp
    src/linalg/matrix.cpp
    src/linalg/vector.cpp
)

target_include_directories(cfd_linalg PUBLIC include)

add_executable(test_gaussian_elimination
    tests/unit_tests/test_gaussian_elimination.cpp
)

target_link_libraries(test_gaussian_elimination PRIVATE cfd_linalg)

add_test(
    NAME test_gaussian_elimination
    COMMAND test_gaussian_elimination
)