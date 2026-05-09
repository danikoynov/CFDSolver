#pragma once

#include "setups/base.hpp"
#include "setups/airfoils/airfoil.hpp"

#include <string>

namespace cfd::setups {

    class AirfoilFlowSetup final : public BaseSetup {
        public:
            AirfoilFlowSetup(
                const std::string& naca_code,
                double chord,
                double freestream_velocity,
                double dx,
                double angle_of_attack_deg
            );

            void impose_boundary_conditions(Grid& grid) const override;

        private:
            bool is_inside_airfoil(
                int i,
                int j,
                int w,
                int h
            ) const;

            void make_solid_cell(
                int i,
                int j,
                Grid& grid,
                BoundaryConditions& bc
            ) const;

        private:
            Airfoil airfoil_;

            double v_air_;
            double dx_;
            double aoa_deg_;
    };

}