#pragma once
#include "linalg/matrix.hpp"
#include "linalg/vector.hpp"
#include "linalg/linear_operator.hpp"

namespace cfd::linalg {

    Vector conjugate_gradient(const Matrix& A, const Vector& b);
}
