#pragma once
#include "setups/base.hpp"

namespace cfd::setups {
    class LidDrivenCavity : public BaseSetup {
        public:
            LidDrivenCavity(double lid_velocity);
            
            void impose_boundary_conditions(
                Grid& grid
            ) const override;

        private:
            double lid_velocity_;
    };

}