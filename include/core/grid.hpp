#pragma once
#include "core/boundary_conditions.hpp"
#include "core/velocity_field.hpp"
#include "core/pressure_field.hpp"

namespace cfd {

    class Grid {

        private:
            std::size_t width_, height_;
            double resolution_;
            VelocityField velocity_;
            PressureField pressure_;
            BoundaryConditions bc_; 
        
        public:
            Grid(std::size_t width, std::size_t height, 
                double resolution, double outside_pressure);

            VelocityField& velocity() { return velocity_; }
            const VelocityField& velocity() const { return velocity_; } 
            
            PressureField& pressure() { return pressure_; }
            const PressureField& pressure() const { return pressure_; }

            BoundaryConditions& boundary_conditions() { return bc_; }
            const BoundaryConditions& boundary_conditions() const { return bc_; }

            std::size_t width() const { return width_; }
            std::size_t height() const { return height_; }
            double resolution() const { return resolution_; }


    };
}