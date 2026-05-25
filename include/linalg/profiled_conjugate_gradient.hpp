#pragma once
#include "linalg/matrix.hpp"
#include "linalg/vector.hpp"
#include "linalg/poisson_operator.hpp"

namespace cfd::linalg {

    Vector profiled_conjugate_gradient_noic(const PoissonOperator& A, const Vector& b);
}
