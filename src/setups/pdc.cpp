#include "pdc.hpp" 


namespace cfd::setups {

    PressureDrivenChannel::PressureDrivenChannel(
        double pressure_difference
    ): pressure_difference_(pressure_difference) {

    }

    void PressureDrivenChannel::impose_boundary_conditions(
        Grid& grid
    ) const {
        std::size_t w = grid.width();
        std::size_t h = grid.height();

        auto& bc = grid.boundary_conditions();

        for (int i = 0; i < static_cast<int>(w); i ++) {
            bc.set_cell_type(i, 0, cfd::CellType::SOLID);
            bc.set_cell_type(i, static_cast<int>(h) - 1, cfd::CellType::SOLID);

            bc.prescribe_u_value(i + 1, 0, 0.0);
            bc.prescribe_u_value(i + 1, static_cast<int>(h) - 1, 0.0);

            bc.prescribe_v_value(i, 1, 0.0);
            bc.prescribe_v_value(i, static_cast<int>(h) - 1, 0.0);
        }

        // moving lid
        for (int i = 2; i < static_cast<int>(w) - 1; i ++) {
            bc.prescribe_u_value(i, static_cast<int>(h) - 2, pressure_difference_);
        }

        for (int j = 0; j < static_cast<int>(h); j ++) {
            bc.set_cell_type(0, j, cfd::CellType::SOLID);
            bc.set_cell_type(static_cast<int>(w) - 1, j, cfd::CellType::SOLID);

            bc.prescribe_v_value(0, j + 1, 0.0);
            bc.prescribe_v_value(static_cast<int>(w) - 1, j + 1, 0.0);

            bc.prescribe_u_value(1, j, 0.0);
            bc.prescribe_u_value(static_cast<int>(w) - 1, j, 0.0);
        }
    }
}