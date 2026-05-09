#include "setups/base.hpp"

namespace cfd::setups {

    void BaseSetup::set_to_solid(
        int i, 
        int j,
        BoundaryConditions& bc
    ) const {
        bc.set_cell_type(i, j, cfd::CellType::SOLID);

        bc.prescribe_u_value(i, j, 0);
        bc.prescribe_u_value(i + 1, j, 0);

        bc.prescribe_v_value(i, j, 0);
        bc.prescribe_v_value(i, j + 1, 0);
    }
}