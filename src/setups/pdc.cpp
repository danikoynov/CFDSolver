#include "setups/pdc.hpp"

namespace cfd::setups {

    PressureDrivenChannel::PressureDrivenChannel(
        double pressure_difference
    )
        : pressure_difference_(pressure_difference) {

    }

    void PressureDrivenChannel::impose_boundary_conditions(
        Grid& grid
    ) const {
        const int w = static_cast<int>(grid.width());
        const int h = static_cast<int>(grid.height());

        auto& bc = grid.boundary_conditions();

        // Corner horizontal velocities.
        bc.prescribe_u_value(0, 0, 0.0);
        bc.prescribe_u_value(0, h - 1, 0.0);

        // Top and bottom solid walls.
        for (int i = 0; i < w; i++) {
            bc.set_cell_type(i, 0, cfd::CellType::SOLID);
            bc.set_cell_type(i, h - 1, cfd::CellType::SOLID);

            // No-slip horizontal velocity on bottom/top walls.
            bc.prescribe_u_value(i + 1, 0, 0.0);
            bc.prescribe_u_value(i + 1, h - 1, 0.0);

            // No-penetration vertical velocity at bottom/top walls.
            bc.prescribe_v_value(i, 1, 0.0);
            bc.prescribe_v_value(i, h - 1, 0.0);
        }

        // Pressure inlet/outlet.
        // Left ghost pressure = pressure_difference_
        // Right ghost pressure = 0
        for (int j = 1; j < h - 1; j++) {
            bc.prescribe_p_value(-1, j, pressure_difference_);
            bc.prescribe_p_value(w, j, 0.0);
        }
    }

}