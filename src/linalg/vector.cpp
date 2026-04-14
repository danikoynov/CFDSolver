#include "linalg/vector.hpp"
#include <stdexcept>

namespace cfd::linalg {

    Vector::Vector(std::size_t n):
        n_(n), data_(n) {
    }

    double& Vector::operator() (std::size_t i){
        return data_[i];
    }

    const double& Vector::operator() (std::size_t i) const{
        return data_[i];
    }

    std::size_t Vector::n() const {
        return n_;
    }

    Vector Vector::operator+(const Vector& other) const {
        if (n_ != other.n_) {
            throw std::invalid_argument("Vector sizes do not match for addition.");
        }

        Vector res(n_);
        for (std::size_t i = 0; i < n_; i++) {
            res(i) = data_[i] + other(i);
        }
        return res;
    }

    Vector Vector::operator-(const Vector& other) const {
        if (n_ != other.n_) {
            throw std::invalid_argument("Vector sizes do not match for subtraction.");
        }

        Vector res(n_);
        for (std::size_t i = 0; i < n_; i++) {
            res(i) = data_[i] - other(i);
        }
        return res;
    }

    Vector& Vector::operator+=(const Vector& other) {
        if (n_ != other.n_) {
            throw std::invalid_argument("Vector sizes do not match for addition.");
        }

        for (std::size_t i = 0; i < n_; i++) {
            data_[i] += other(i);
        }
        return *this;
    }

    Vector& Vector::operator-=(const Vector& other) {
        if (n_ != other.n_) {
            throw std::invalid_argument("Vector sizes do not match for subtraction.");
        }

        for (std::size_t i = 0; i < n_; i++) {
            data_[i] -= other(i);
        }
        return *this;
    }



    double Vector::dot(const Vector& a, const Vector& b) {
        if (a.n_ != b.n_) {
            throw std::invalid_argument("Vector sizes do not match for dot product.");
        }

        double res = 0.0;
        for (std::size_t i = 0; i < a.n_; i++) {
            res += a.data_[i] * b.data_[i];
        }
        return res;
    }


    Vector Vector::operator*(double alpha) const {
        Vector res(n_);
        for (std::size_t i = 0; i < n_; i++) {
            res(i) = data_[i] * alpha;
        }
        return res;
    }

    Vector operator*(double alpha, const Vector& v) {
        return v * alpha;
    }

}