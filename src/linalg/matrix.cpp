#pragma once
#include "linalg/matrix.hpp"
#include <utility>


namespace cfd::linalg {

    Matrix::Matrix(std::size_t rows, std::size_t cols):
        rows_(rows), cols_(cols), data_(rows * cols) {
    }

    double& Matrix::operator() (std::size_t i, std::size_t j){
        return data_[i * cols_ + j];
    }

    const double& Matrix::operator() (std::size_t i, std::size_t j) const{
        return data_[i * cols_ + j];
    }

    void Matrix::swap_rows(std::size_t i, std::size_t j) {
        for (std::size_t k = 0; k < cols_; k ++) {
            std::swap((*this)(i, k), (*this)(j, k));
        }
    }

    std::size_t Matrix::rows() const {
        return rows_;
    }

    std::size_t Matrix::cols() const {
        return cols_;
    }
}