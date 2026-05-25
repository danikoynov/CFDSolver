#pragma once

#include "linalg/matrix.hpp"
#include "linalg/vector.hpp"
#include "linalg/poisson_operator.hpp"

#include <filesystem>
#include <optional>

namespace cfd::linalg {

    Vector profiled_conjugate_gradient(
        const PoissonOperator& A,
        const Vector& b,
        bool save_data,
        std::optional<std::filesystem::path> output_dir
    );

}