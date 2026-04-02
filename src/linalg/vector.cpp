#include "linalg/vector.hpp"


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
}