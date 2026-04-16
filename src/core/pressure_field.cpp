#include "core/pressure_field.hpp"

#include <stdexcept>

namespace cfd {

    PressureField::PressureField(
        std::size_t width,
        std::size_t height
    )
        : width_(width),
        height_(height),
        p_(width * height) {}

    void PressureField::check_bounds(int i, int j) const {
        // Throws if (i, j) is outside the pressure field bounds
        
        const int width = static_cast<int>(width_);
        const int height = static_cast<int>(height_);

        if (i < 0 || i >= width ||
            j < 0 || j >= height) {
            throw std::out_of_range("pressure indexes out of bounds");
        }
    }

    double& PressureField::get_p(int i, int j) {
        // Returns access to the pressure at cell [i, j]

        check_bounds(i, j);

        return p_[static_cast<std::size_t>(i) * height_ +
                static_cast<std::size_t>(j)];
    }

    const double& PressureField::get_p(int i, int j) const {
        // Returns read-only access to the pressure at cell [i, j]

        check_bounds(i, j);

        std::size_t id = static_cast<std::size_t>(i) * height_ +
                static_cast<std::size_t>(j);
                
        return p_[id];
    }

}  