#pragma once
#include "linalg/matrix.hpp"
#include "linalg/vector.hpp"
#include "linalg/linear_operator.hpp"

namespace cfd::linalg {

    Vector conjugate_gradient_profiled(const LinearOperator& A, const Vector& b);
}
