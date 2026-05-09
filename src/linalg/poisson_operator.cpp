#include "linalg/poisson_operator.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfd::linalg {

    PoissonOperator::PoissonOperator(
        std::size_t width,
        std::size_t height
    ) :
        width_(width),
        height_(height),
        size_(width * height),
        assembly_rows_(size_),
        row_ptr_(size_ + 1) {

        if (width_ == 0 || height_ == 0) {
            throw std::runtime_error("PoissonOperator requires nonzero width and height.");
        }
    }

    Vector PoissonOperator::apply(const Vector& v) const {
        if (v.n() != size_) {
            throw std::runtime_error("Input vector size does not match Poisson operator size.");
        }

        Vector res(size_);

        for (std::size_t i = 0; i < size_; i++) {
            for (std::size_t j = row_ptr_[i]; j < row_ptr_[i + 1]; j++) {
                const auto& [col, coeff] = csr_entries_[j];
                res(i) += coeff * v(col);
            }
        }

        return res;
    }

    double& PoissonOperator::operator()(int x, int y) {
        if (x < 0 || x >= static_cast<int>(size_)) {
            throw std::range_error("Row index outside of Poisson operator range.");
        }

        if (y < 0 || y >= static_cast<int>(size_)) {
            throw std::range_error("Column index outside of Poisson operator range.");
        }

        std::size_t row = static_cast<std::size_t>(x);
        std::size_t col = static_cast<std::size_t>(y);

        std::size_t pivot = 0;
        while (
            pivot < assembly_rows_[row].size() &&
            assembly_rows_[row][pivot].first != col
        ) {
            pivot++;
        }

        if (pivot == assembly_rows_[row].size()) {
            assembly_rows_[row].emplace_back(col, 0.0);
        }

        return assembly_rows_[row][pivot].second;
    }

    const double& PoissonOperator::operator()(int x, int y) const {
        if (x < 0 || x >= static_cast<int>(size_)) {
            throw std::range_error("Row index outside of Poisson operator range.");
        }

        if (y < 0 || y >= static_cast<int>(size_)) {
            throw std::range_error("Column index outside of Poisson operator range.");
        }

        std::size_t row = static_cast<std::size_t>(x);
        std::size_t col = static_cast<std::size_t>(y);

        for (const auto& [entry_col, coeff] : assembly_rows_[row]) {
            if (entry_col == col) {
                return coeff;
            }
        }

        throw std::runtime_error("Accessed Poisson operator entry does not exist.");
    }

    double PoissonOperator::get(int x, int y) const {
        if (x < 0 || x >= static_cast<int>(size_)) {
            return 0.0;
        }

        if (y < 0 || y >= static_cast<int>(size_)) {
            return 0.0;
        }

        std::size_t row = static_cast<std::size_t>(x);
        std::size_t col = static_cast<std::size_t>(y);

        auto it = std::find_if(
            assembly_rows_[row].begin(),
            assembly_rows_[row].end(),
            [col](const auto& entry) {
                return entry.first == col;
            }
        );

        if (it == assembly_rows_[row].end()) {
            return 0.0;
        }

        return it->second;
    }

    double PoissonOperator::get_diag(int k) const {
        return get(k, k);
    }

    double PoissonOperator::get_bottom(int k) const {
        return get(k, k - static_cast<int>(width_));
    }

    double PoissonOperator::get_top(int k) const {
        return get(k, k + static_cast<int>(width_));
    }

    double PoissonOperator::get_left(int k) const {
        if (k < 0 || k >= static_cast<int>(size_)) {
            return 0.0;
        }

        if (k % static_cast<int>(width_) == 0) {
            return 0.0;
        }

        return get(k, k - 1);
    }

    double PoissonOperator::get_right(int k) const {
        if (k < 0 || k >= static_cast<int>(size_)) {
            return 0.0;
        }

        if ((k + 1) % static_cast<int>(width_) == 0) {
            return 0.0;
        }

        return get(k, k + 1);
    }

    void PoissonOperator::sterilize() {
        for (auto& row : assembly_rows_) {
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

    void PoissonOperator::init_csr() {
        csr_entries_.clear();

        for (std::size_t i = 0; i < size_; i++) {
            row_ptr_[i] = csr_entries_.size();

            for (const auto& [col, coeff] : assembly_rows_[i]) {
                csr_entries_.emplace_back(col, coeff);
            }
        }

        row_ptr_[size_] = csr_entries_.size();
    }

    void PoissonOperator::init_preconditioner() {
        ic_inv_diag_.assign(size_, 0.0);
        ic_bottom_.assign(size_, 0.0);
        ic_left_.assign(size_, 0.0);

        constexpr double SAFETY = 1e-12;

        const int n = static_cast<int>(size_);
        const int w = static_cast<int>(width_);

        for (int k = 0; k < n; k++) {
            double diag = get(k, k);
            double e = diag;

            int bottom = k - w;
            if (bottom >= 0) {
                double a_bottom = get(k, bottom);
                double l_bottom = a_bottom * ic_inv_diag_[bottom];

                ic_bottom_[k] = l_bottom;
                e -= l_bottom * l_bottom;
            }

            int left = k - 1;
            if (k % w != 0 && left >= 0) {
                double a_left = get(k, left);
                double l_left = a_left * ic_inv_diag_[left];

                ic_left_[k] = l_left;
                e -= l_left * l_left;
            }

            if (e <= SAFETY) {
                throw std::runtime_error(
                    "Incomplete Cholesky breakdown: non-positive diagonal."
                );
            }

            ic_inv_diag_[k] = 1.0 / std::sqrt(e);
        }
    }

    Vector PoissonOperator::apply_preconditioner(const Vector& r) const {
        if (r.n() != size_) {
            throw std::runtime_error("Residual size does not match Poisson operator size.");
        }

        if (
            ic_inv_diag_.size() != size_ ||
            ic_bottom_.size() != size_ ||
            ic_left_.size() != size_
        ) {
            throw std::runtime_error("Preconditioner has not been initialized.");
        }

        Vector q(size_);
        Vector z(size_);

        const int n = static_cast<int>(size_);
        const int w = static_cast<int>(width_);

        // Forward solve: L q = r
        for (int k = 0; k < n; k++) {
            double t = r(k);

            int bottom = k - w;
            if (bottom >= 0) {
                t -= ic_bottom_[k] * q(bottom);
            }

            int left = k - 1;
            if (k % w != 0 && left >= 0) {
                t -= ic_left_[k] * q(left);
            }

            q(k) = t * ic_inv_diag_[k];
        }

        // Backward solve: L^T z = q
        for (int k = n - 1; k >= 0; k--) {
            double t = q(k);

            int top = k + w;
            if (top < n) {
                t -= ic_bottom_[top] * z(top);
            }

            int right = k + 1;
            if ((k + 1) % w != 0 && right < n) {
                t -= ic_left_[right] * z(right);
            }

            z(k) = t * ic_inv_diag_[k];
        }

        return z;
    }

    void PoissonOperator::finalize() {
        sterilize();
        init_csr();
        init_preconditioner();
    }

}