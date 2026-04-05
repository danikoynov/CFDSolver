#include "core/pressure_field.hpp"
#include <stdexcept>

namespace cfd {

    PressureField::PressureField(
        std::size_t width, std::size_t height, 
        double outside_pressure
    ): width_(width), height_(height),
        outside_pressure_(outside_pressure), p_(width * height) {}

    void PressureField::check_bounds(int i, int j) const
    {
        if (i < 0 || i >= static_cast<int>(width_) ||
            j < 0 || j >= static_cast<int>(height_)) {
            throw std::out_of_range("pressure indexes out of bounds");
        }
    }

    double& PressureField::get_p(int i, int j) {
        check_bounds(i, j);
        return p_[i * height_ + j];
    }

    const double& PressureField::get_p(int i, int j) const {
        check_bounds(i, j);
        return p_[i * height_ + j];
    }

    double PressureField::read_p_or_outside(int i, int j) const {
        if (i < 0 || i >= static_cast<int>(width_) ||
            j < 0 || j >= static_cast<int>(height_)) 
            return outside_pressure_;
        
        return get_p(i, j);
    }

    double PressureField::outside_pressure() const {
        return outside_pressure_;
    }


}