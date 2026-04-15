#include "core/boundary_conditions.hpp"
#include <stdexcept>

namespace cfd {

    BoundaryConditions::BoundaryConditions(std::size_t width, std::size_t height)
    : width_(width),
      height_(height),
      type_(width * height, CellType::FLUID) {
    }

    void BoundaryConditions::prescribe_u_value(int i, int j, double val) {
        /// Prescribe a value for u in cell [i, j]
        
        const int width = static_cast<int>(width_);
        const int height = static_cast<int>(height_);

        if (i < 0 || i > width ||
            j < 0 || j >= height) {
            throw std::out_of_range("boundary condition u indexes out of bounds");
        }

        const std::size_t id =
            static_cast<std::size_t>(i) * height_ +
            static_cast<std::size_t>(j);

        const bool key_present =
            prescribed_u_.find(id) != prescribed_u_.end();

        if (key_present && prescribed_u_.at(id) != val) {
            throw std::runtime_error(
                "boundary condition for u already has a different value prescribed"
            );
        }

        if (!key_present) {
            prescribed_u_.insert({id, val});
        }
    }

    void BoundaryConditions::prescribe_v_value(int i, int j, double val) {
        /// Prescribe a value for v in cell [i, j]

        const int width = static_cast<int>(width_);
        const int height = static_cast<int>(height_);

        if (i < 0 || i >= width ||
            j < 0 || j > height) {
            throw std::out_of_range("boundary condition v indexes out of bounds");
        }

        const std::size_t id =
            static_cast<std::size_t>(i) * (height_ + 1) +
            static_cast<std::size_t>(j);

        const bool key_present =
            prescribed_v_.find(id) != prescribed_v_.end();

        if (key_present && prescribed_v_.at(id) != val) {
            throw std::runtime_error(
                "boundary condition for v already has a different value prescribed"
            );
        }

        if (!key_present) {
            prescribed_v_.insert({id, val});
        }
    }

    void BoundaryConditions::prescribe_p_value(int i, int j, double val) {
        /// Prescribe a value for p in cell [i, j]
        /// The cell could be on the boundary of the grid

        const int width = static_cast<int>(width_);
        const int height = static_cast<int>(height_);

        if (i < -1 || i > width ||
            j < -1 || j > height_) {
            throw std::out_of_range("boundary condition p indexes out of bounds");
        }

        const int id = i * static_cast<int>(height_) + j;

        const bool key_present =
            prescribed_p_.find(id) != prescribed_p_.end();

        if (key_present && prescribed_p_.at(id) != val) {
            throw std::runtime_error(
                "boundary condition for p already has a different value prescribed"
            );
        }

        if (!key_present) {
            prescribed_p_.insert({id, val});
        }
    }


    const std::unordered_map<std::size_t, double>& BoundaryConditions::prescribed_u() const {
        return prescribed_u_;
    }

    const std::unordered_map<std::size_t, double>& BoundaryConditions::prescribed_v() const {
        return prescribed_v_;
    }
    
    const std::unordered_map<int, double>& BoundaryConditions::prescribed_p() const {
        return prescribed_p_;
    }
    
    double BoundaryConditions::prescribed_u(int i, int j) const {
        /// Returns prescribed u value on left face of cell [i, j], 
        /// throws if it is not prescribed
                            
        const int width = static_cast<int>(width_);
        const int height = static_cast<int>(height_);

        if (i < 0 || i > width ||
            j < 0 || j >= height) {
            throw std::out_of_range("boundary condition u indexes out of bounds");
        }

        const std::size_t id =
            static_cast<std::size_t>(i) * height_ +
            static_cast<std::size_t>(j);

        return prescribed_u_.at(id);
    }

    double BoundaryConditions::prescribed_v(int i, int j) const {
        /// Returns prescribed v value on down face of cell [i, j], 
        /// throws if it is not prescribed
                    
        const int width = static_cast<int>(width_);
        const int height = static_cast<int>(height_);

        if (i < 0 || i >= width ||
            j < 0 || j > height) {
            throw std::out_of_range("boundary condition v indexes out of bounds");
        }

        const std::size_t id =
            static_cast<std::size_t>(i) * (height_ + 1) +
            static_cast<std::size_t>(j);

        return prescribed_v_.at(id);
    }

    double BoundaryConditions::prescribed_p(int i, int j) const {
        /// Returns prescribed p value in cell [i, j], 
        /// throws if it is not prescribed
                
        const int width = static_cast<int>(width_);
        const int height = static_cast<int>(height_);

        if (i < -1 || i > width ||
            j < -1 || j > height) {
            throw std::out_of_range("boundary condition p indexes out of bounds");
        }

        const int id = i * height + j;

        const bool key_present =
            prescribed_p_.find(id) != prescribed_p_.end();

        if (!key_present) {
            return DEFAULT_PRESSURE;
        }

        return prescribed_p_.at(id);
    }

    bool BoundaryConditions::is_u_prescribed(int i, int j) const {
        /// Returns whether or not left face of cell [i, j] has prescribed u value

        const int width = static_cast<int>(width_);
        const int height = static_cast<int>(height_);

        if (i < 0 || i > width ||
            j < 0 || j >= height) {
            throw std::out_of_range("boundary condition u indexes out of bounds");
        }

        const std::size_t id =
            static_cast<std::size_t>(i) * height_ +
            static_cast<std::size_t>(j);

        return prescribed_u_.find(id) != prescribed_u_.end();
    }

    bool BoundaryConditions::is_v_prescribed(int i, int j) const {
        /// Returns whether or not down face of cell [i, j] has prescribed v value
        
        const int width = static_cast<int>(width_);
        const int height = static_cast<int>(height_);

        if (i < 0 || i >= width ||
            j < 0 || j > height) {
            throw std::out_of_range("boundary condition v indexes out of bounds");
        }

        const std::size_t id =
            static_cast<std::size_t>(i) * (height_ + 1) +
            static_cast<std::size_t>(j);

        return prescribed_v_.find(id) != prescribed_v_.end();
    }

    bool BoundaryConditions::is_p_prescribed(int i, int j) const {
        /// Returns whether or not cell [i, j] has prescribed p value
        
        const int width = static_cast<int>(width_);
        const int height = static_cast<int>(height_);

        if (i < 0 || i >= width ||
            j < 0 || j >= height) {
            throw std::out_of_range("boundary condition p indexes out of bounds");
        }

        const std::size_t id =
            static_cast<std::size_t>(i) * height_ +
            static_cast<std::size_t>(j);

        return prescribed_p_.find(id) != prescribed_p_.end();
    }

    CellType BoundaryConditions::type(int i, int j) const {
        /// Return the type of cell [i, j]

        const int width = static_cast<int>(width_);
        const int height = static_cast<int>(height_);

        if (i < 0 || i >= width ||
            j < 0 || j >= height) {
            return CellType::BOUNDARY;
        }

        const std::size_t id =
            static_cast<std::size_t>(i) * height_ +
            static_cast<std::size_t>(j);

        return type_[id];
    }

    void BoundaryConditions::set_cell_type(int i, int j, CellType cell_type) {
        /// Set the type of cell [i, j]
        
        const int width = static_cast<int>(width_);
        const int height = static_cast<int>(height_);

        if (i < 0 || i >= width ||
            j < 0 || j >= height) {
            throw std::out_of_range("setting cell type indexes out of bounds");
        }

        const std::size_t id =
            static_cast<std::size_t>(i) * height_ +
            static_cast<std::size_t>(j);

        type_[id] = cell_type;
    }


}
