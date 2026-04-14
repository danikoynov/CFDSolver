#pragma once
#include <optional>
#include "core/grid.hpp"
#include "linalg/gaussian_elimination.hpp"
#include "linalg/linear_operator.hpp"
#include "linalg/matrix.hpp"
#include "linalg/vector.hpp"

/**
 * Workflow of Simulator tick
 *  1. Determine timestep
 *  2. Advect Velcoity Field
 *  3. Apply Boundary Conditions
 *  4. Apply Body Forces
 *  5. Apply Boundary Conditions
 *  6. Project Velocity Field using Pressure field
 *  7. Apply Boundary Conditions
 */

namespace cfd {
    
    class Simulator {
        
        static constexpr double FALLBACK_MAX_TIMESTEP = 0.1;
        static constexpr double CFL = 0.5;
        static constexpr double TOL = 1e-9; /// General Numeric Tolerance
        static constexpr double GRAVITY_AY = -9.80665; 

        private:
            Grid grid_;
            double fluid_density_;
            double viscosity_;
            bool apply_gravity_;

        public:
            Simulator(std::size_t width, std::size_t height, 
                double resolution, double outside_pressure, 
                double fluid_density, bool apply_gravity,
                double viscosity);

            double determine_timestep() const;
            void advect(double timestep);
            void apply_body_forces(double timestpe, double ax, double ay);
            void apply_boundary_conditions();
            void project(double timestep);
            void apply_viscosity(double timestep);
            std::optional<double> build_equation(int i, int j, double timestep,
                linalg::LinearOperator& poisson_matrix, linalg::Vector& poisson_rhs);
            void apply_pressure_gradient(double timestep);
            void set_pressure_values(const linalg::Vector& pressure_values);
            std::optional<double> get_boundary_pressure_constraint(int i, int j, double timestep);
            void tick();
            const Grid& grid() const;
            Grid& grid();
    };
}