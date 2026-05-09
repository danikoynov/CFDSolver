#pragma once
#include "core/grid.hpp"

namespace cfd::setups {

    class BaseSetup {
        public:
            BaseSetup() {}

            virtual void impose_boundary_conditions(
                Grid& grid
            ) const = 0;
    
            void set_to_solid(
                int i,
                int j,
                BoundaryConditions& bc
            ) const;

    };

}