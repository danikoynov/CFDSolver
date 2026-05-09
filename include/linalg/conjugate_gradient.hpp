#pragma once
#include "linalg/matrix.hpp"
#include "linalg/vector.hpp"
#include "linalg/poisson_operator.hpp"

namespace cfd::linalg {

    Vector conjugate_gradient(const PoissonOperator& A, const Vector& b);
}
