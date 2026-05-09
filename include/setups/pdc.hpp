#pragma once
#include "base.hpp"

namespace cfd::setups {

    class PressureDrivenChannel : public BaseSetup {
        public:
            PressureDrivenChannel(double pressure_difference_);

            void impose_boundary_conditions(
                Grid& grid
            ) const override;

        private:
            double pressure_difference_;
            
    };
}