#pragma once
#include "vector.hpp"
#include "matrix.hpp"

namespace cfd::linalg {

    Vector gaussian_elimination(Matrix A, Vector b);
}