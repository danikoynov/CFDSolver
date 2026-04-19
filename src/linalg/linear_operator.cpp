#include "linalg/linear_operator.hpp"
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace cfd::linalg {

    LinearOperator::LinearOperator(
        std::size_t input_size, std::size_t output_size) :
        input_size_(input_size), 
        output_size_(output_size),
        dep_(output_size),
        row_ptr_(output_size + 1) {
    }

    void LinearOperator::init_csr() {
        csr_dep_.clear();

        for (std::size_t i = 0; i < output_size_; i++) {
            row_ptr_[i] = csr_dep_.size();
            for (const auto& [col, coeff] : dep_[i]) {
                csr_dep_.emplace_back(col, coeff);
            }
        }

        row_ptr_[output_size_] = csr_dep_.size();
    }

    Vector LinearOperator::apply(const Vector& v) const {

        if (input_size_ != v.n()) {
            throw std::runtime_error("Input vector size does not match linear operator input size."); 
        }
        
        Vector res(output_size_);

        for (std::size_t i = 0; i < output_size_; i ++) {
            for (std::size_t j = row_ptr_[i]; j < row_ptr_[i + 1]; j ++) {
                const auto& [col, coeff] = csr_dep_[j];
                res(i) += coeff * v(col);
            }
        }

        return res;
    }

   double& LinearOperator::operator()(int x, int y) {
        if (x < 0 || x >= static_cast<int>(output_size_)) {
            throw std::range_error("Dependent variable index outside of range");
        }
        if (y < 0 || y >= static_cast<int>(input_size_)) {
            throw std::range_error("Independent variable index outside of range");
        }

        std::size_t row = static_cast<std::size_t>(x);
        std::size_t col = static_cast<std::size_t>(y);

        std::size_t pivot = 0;
        while (pivot < dep_[row].size() && dep_[row][pivot].first != col) {
            pivot++;
        }

        if (pivot == dep_[row].size()) {
            dep_[row].emplace_back(col, 0.0);
        }

        return dep_[row][pivot].second;
    }

    const double& LinearOperator::operator()(int x, int y) const { 
        if (x < 0 || x >= output_size_) {
            throw std::range_error("Dependent variable index outside of range");
        }
        if (y < 0 || y >= input_size_) {
            throw std::range_error("Independent variable index outside of range");
        }

        std::size_t pivot = 0;
        while(pivot < dep_[x].size() && dep_[x][pivot].first != y) {
            pivot ++;
        }

        if (pivot != dep_[x].size())
            return dep_[x][pivot].second;

        throw std::runtime_error("Accessed dependency does not exist");
    }

    void LinearOperator::sterilize() {
        for (auto& row : dep_) {
            row.erase(
                std::remove_if(
                    row.begin(),
                    row.end(),
                    [](const auto& entry) {
                        return std::abs(entry.second) < EPS;
                    }
                ),
                row.end()
            );
        }
    }
}