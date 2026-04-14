#include "linalg/linear_operator.hpp"
#include <cmath>
#include <stdexcept>

namespace cfd::linalg {

    LinearOperator::LinearOperator(
        const Matrix& A) :
        input_size_(A.cols()), output_size_(A.rows()) {

        std::size_t r = A.rows();
        std::size_t c = A.cols();
        
        dep_.resize(output_size_);

        for (std::size_t i = 0; i < r; i ++) {
            for (std::size_t j = 0; j < c; j ++) {
                
                if (std::abs(A(i, j)) < EPS) {
                    continue;
                }

                dep_[i].push_back(std::make_pair(j, A(i, j)));
            }
        }
    }

    Vector LinearOperator::apply(const Vector& v) const {

        if (input_size_ != v.n()) {
            throw std::runtime_error("Input vector size does not match linear operator input size."); 
        }
        
        Vector res(output_size_);

        for (std::size_t i = 0; i < output_size_; i ++) {
            for (const auto& [col, coeff] : dep_[i]) {
                res(i) += coeff * v(col);
            }
        }

        return res;
    }

}