#include "core/boundary_conditions.hpp"
#include <stdexcept>

namespace cfd {

    BoundaryConditions::BoundaryConditions(
        std::size_t width, std::size_t height
    ):  width_(width), height_(height), type_(width * height, CellType::FLUID) {}


    void BoundaryConditions::prescribe_u_value(int i, int j, double val) {
        
        if (i < 0 || i >  static_cast<int>(width_) ||
            j < 0 || j >= static_cast<int>(height_)) {
            throw std::out_of_range("boundary condition u indexes out of bounds");
        }

        std::size_t id = 
            static_cast<std::size_t>(i) * height_ + 
            static_cast<std::size_t>(j);

        if (prescribed_u_.find(id) != prescribed_u_.end()) {
            throw std::runtime_error("boundary condition for u already prescribed");
        }

        prescribed_u_.insert({id, val});
    }

    void BoundaryConditions::prescribe_v_value(int i, int j, double val) {
        
        if (i < 0 || i >=  static_cast<int>(width_) ||
            j < 0 || j >   static_cast<int>(height_)) {
            throw std::out_of_range("boundary condition v indexes out of bounds");
        }

        std::size_t id = 
            static_cast<std::size_t>(i) * (height_ + 1) + 
            static_cast<std::size_t>(j);

        if (prescribed_v_.find(id) != prescribed_v_.end()) {
            throw std::runtime_error("boundary condition for v already prescribed");
        }

        prescribed_v_.insert({id, val});
    }


    const std::unordered_map<std::size_t, double>& BoundaryConditions::prescribed_u() const {
        return prescribed_u_;
    }

    const std::unordered_map<std::size_t, double>& BoundaryConditions::prescribed_v() const {
        return prescribed_v_;
    }
}