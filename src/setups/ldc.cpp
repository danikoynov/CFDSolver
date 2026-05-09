#include "setups/ldc.hpp"

namespace cfd::setups {

    LidDrivenCavity::LidDrivenCavity(
        double lid_velocity
    ):
        lid_velocity_(lid_velocity) {

    }

    void LidDrivenCavity::impose_boundary_conditions(
        Grid& grid
    ) const {
        
        std::size_t w = grid.width();
        std::size_t h = grid.height();

        auto& bc = grid.boundary_conditions();
        for (int i = 0; i < w; i ++) {
            set_to_solid(i, 0, bc);
            set_to_solid(i, h - 1, bc);
        }
        
        for (int j = 0; j < h; j ++) {
            set_to_solid(0, j, bc);
            set_to_solid(w - 1, j, bc);
        }

        for (int i = 2; i < w - 1; i ++) {
            bc.prescribe_u_value(i, h - 2, lid_velocity_);
            bc.prescribe_u_value(i, 1, 0);
        }

        for (int j = 2; j < h - 1; j ++) {
            bc.prescribe_v_value(1, j, 0);
            bc.prescribe_v_value(w - 2, j, 0);
        }
    }   
}