#pragma once
#include "linalg/matrix.hpp"
#include "linalg/vector.hpp"
#include <vector>

namespace cfd::linalg {
    
    inline constexpr double EPS = 1e-9;

    class LinearOperator {

        private:
            std::size_t input_size_;
            std::size_t output_size_;
            std::vector<std::vector<std::pair<std::size_t, double>>> dep_;

        public:
        
            LinearOperator(const Matrix& A);
            Vector apply(const Vector& v) const;
    };

}