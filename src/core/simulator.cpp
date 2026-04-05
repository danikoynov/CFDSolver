#include "core/simulator.hpp"
#include "linalg/gaussian_elimination.hpp"
#include <optional>
#include <stdexcept>
#include <iostream> /// remove

namespace cfd {

    Simulator::Simulator(std::size_t width, std::size_t height,
        double resolution, double outside_pressure, double fluid_density,
        bool apply_gravity):
        grid_(width, height, resolution, outside_pressure), 
        fluid_density_(fluid_density), apply_gravity_(apply_gravity) {}

    double Simulator::determine_timestep() const {
        double max_vel = grid_.velocity().get_max_velocity_component();

        if (max_vel == 0.0)
            return FALLBACK_MAX_TIMESTEP;
    
        return std::min(FALLBACK_MAX_TIMESTEP, CFL * grid_.resolution() / max_vel);
    }

    void Simulator::advect(double timestep) {
        /// Performs Semi-Lagrangian advection on the velocity field
        const auto& old_vel = grid_.velocity();
        const std::size_t w = grid_.width();
        const std::size_t h = grid_.height();
        const double dx = grid_.resolution();

        VelocityField next_field(w, h, dx);

        /// Advect vertical faces
        for (std::size_t i = 0; i <= w; i ++)
            for (std::size_t j = 0; j < h; j ++) {
                auto velocity = old_vel.sample_at_vertical_face(i, j);

                double cur_x = static_cast<double>(i) * dx;
                double cur_y = (static_cast<double>(j) + 0.5) * dx;
                
                double prev_x = cur_x - timestep * velocity.u;
                double prev_y = cur_y - timestep * velocity.v;

                auto new_velocity = old_vel.sample_at_coordinates(prev_x, prev_y);

                next_field.get_u(i, j) = new_velocity.u;
            }

        /// Advect horizontal faces
        for (std::size_t i = 0; i < w; i ++)
            for (std::size_t j = 0; j <= h; j ++) {
                auto velocity = old_vel.sample_at_horizontal_face(i, j);

                double cur_x = (static_cast<double>(i) + 0.5) * dx;
                double cur_y = (static_cast<double>(j)) * dx;
                
                double prev_x = cur_x - timestep * velocity.u;
                double prev_y = cur_y - timestep * velocity.v;

                auto new_velocity = old_vel.sample_at_coordinates(prev_x, prev_y);

                next_field.get_v(i, j) = new_velocity.v;
            }    
        
        grid_.velocity() = std::move(next_field);
    }

    void Simulator::apply_body_forces(double timestep, double ax, double ay) {
        auto& velocity = grid_.velocity();
        const std::size_t w = grid_.width();
        const std::size_t h = grid_.height();

        /// Apply on vertical faces
        for (std::size_t i = 0; i <= w; i ++)
            for (std::size_t j = 0; j < h; j ++) {
                velocity.get_u(i, j) += timestep * ax;
            }

        /// Apply on horizontal faces
        for (std::size_t i = 0; i < w; i ++)
            for (std::size_t j = 0; j <= h; j ++) {
                velocity.get_v(i, j) += timestep * ay;
            }
    }

    void Simulator::apply_boundary_conditions() {
        const auto& bc = grid_.boundary_conditions();

        const auto& prescribed_u = bc.prescribed_u();
        const auto& prescribed_v = bc.prescribed_v();

        for (const auto& [key, value] : prescribed_u) {
            std::size_t i = key / (grid_.height());
            std::size_t j = key % (grid_.height());
            
            grid_.velocity().get_u(i, j) = value;
        }

        for (const auto& [key, value] : prescribed_v) {
            std::size_t i = key / (grid_.height() + 1);
            std::size_t j = key % (grid_.height() + 1);
            
            grid_.velocity().get_v(i, j) = value;
        }
    }

    std::optional<double> Simulator::get_boundary_pressure_constraint(int i, int j, double timestep) {
        const auto& bc = grid_.boundary_conditions();
        const auto& vf = grid_.velocity();
        const auto& pf = grid_.pressure();
        
        std::optional<double> constr;
        
        auto add_pressure_constraint = [&](double p) {
            if (!constr.has_value()) {
                constr = p;
                return;
            }

            if (std::abs(*constr - p) > TOL) {
                throw std::runtime_error(
                    "Incompatible pressure/boundary conditions:\n"
                    "cell pressure is overconstrained"
                );
            }
        };
        
        /**std::cout<<"Cell: " << i << " | " << j << " LEFT: ";
        if (bc.type(i - 1, j) == CellType::BOUNDARY)
            std::cout << "Boundary" << std::endl;
        else if (bc.type(i - 1, j) == CellType::SOLID)
            std::cout << "Solid" << std::endl;
        else if (bc.type(i - 1, j) == CellType::FLUID)
            std::cout << "Fluid" << std::endl;*/

        const double factor = fluid_density_ * grid_.resolution() / timestep; 
        if (bc.type(i - 1, j) == CellType::BOUNDARY && bc.is_u_prescribed(i, j)) {
            double cell_pressure = pf.outside_pressure() + 
                factor * (vf.get_u(i, j) - bc.prescribed_u(i, j));
            add_pressure_constraint(cell_pressure);
        }

        if (bc.type(i, j - 1) == CellType::BOUNDARY && bc.is_v_prescribed(i, j)) {
            double cell_pressure = pf.outside_pressure() + 
                factor * (vf.get_v(i, j) - bc.prescribed_v(i, j));
            add_pressure_constraint(cell_pressure);
        }

        if (bc.type(i + 1, j) == CellType::BOUNDARY && bc.is_u_prescribed(i + 1, j)) {
            double cell_pressure = pf.outside_pressure() - 
                factor * (vf.get_u(i + 1, j) - bc.prescribed_u(i + 1, j));
            add_pressure_constraint(cell_pressure);
        }

        if (bc.type(i, j + 1) == CellType::BOUNDARY && bc.is_v_prescribed(i, j + 1)) {
            double cell_pressure = pf.outside_pressure() - 
                factor * (vf.get_v(i, j + 1) - bc.prescribed_v(i, j + 1));
            add_pressure_constraint(cell_pressure);
        }

        return constr;
    }

    void Simulator::build_equation(int i, int j, double timestep,
                linalg::Matrix& poisson_matrix, linalg::Vector& poisson_rhs) {
        
        std::size_t h = grid_.height();
        auto convert_index = [h](std::size_t i, std::size_t j) {
            return i * h + j;
        };

        const auto& bc = grid_.boundary_conditions();
        const auto& vf = grid_.velocity();
        const auto& pf = grid_.pressure();
        
        std::size_t cell_id = convert_index(i, j);

        if (bc.type(i, j) == CellType::SOLID) {
            /// We enforce 0 pressure on the solid cells
            poisson_matrix(cell_id, cell_id) = 1;
            poisson_rhs(cell_id) = 0;
            return;
        }
  
        /// Handle edge case with pressure and velocity fixed boundary
        std::optional<double> boundary_constr = get_boundary_pressure_constraint(i , j, timestep);
        
        if (boundary_constr.has_value()) {
            poisson_matrix(cell_id, cell_id) = 1;
            poisson_rhs(cell_id) = *boundary_constr;
            return;
        }

        double div = vf.get_divergence(i, j);

        poisson_rhs(cell_id) = 
            - (fluid_density_ * grid_.resolution() * grid_.resolution() * div) / timestep; // baseline RHS
        
        const double factor = (fluid_density_ * grid_.resolution() / timestep);
        if (bc.type(i - 1, j) == CellType::FLUID) {
            poisson_matrix(cell_id, cell_id) ++;
            poisson_matrix(cell_id, convert_index(i - 1, j)) --;
        }
        else if (bc.type(i - 1, j) == CellType::BOUNDARY) {
            poisson_matrix(cell_id, cell_id) ++;
            poisson_rhs(cell_id) += pf.outside_pressure();
        }
        else if (bc.type(i - 1, j) == CellType::SOLID) {
            double value = - factor * (vf.get_u(i, j) - bc.prescribed_u(i, j));
            poisson_rhs(cell_id) += value;
        }

        if (bc.type(i + 1, j) == CellType::FLUID) {
            poisson_matrix(cell_id, cell_id) ++;
            poisson_matrix(cell_id, convert_index(i + 1, j)) --;
        }
        else if (bc.type(i + 1, j) == CellType::BOUNDARY) {
            poisson_matrix(cell_id, cell_id) ++;
            poisson_rhs(cell_id) += pf.outside_pressure();
        }
        else if (bc.type(i + 1, j) == CellType::SOLID) {
            double value = factor * (vf.get_u(i + 1, j) - bc.prescribed_u(i + 1, j));
            poisson_rhs(cell_id) += value;
        }

        if (bc.type(i, j - 1) == CellType::FLUID) {
            poisson_matrix(cell_id, cell_id) ++;
            poisson_matrix(cell_id, convert_index(i, j - 1)) --;
        }
        else if (bc.type(i, j - 1) == CellType::BOUNDARY) {
            poisson_matrix(cell_id, cell_id) ++;
            poisson_rhs(cell_id) += pf.outside_pressure();
        }
        else if (bc.type(i, j - 1) == CellType::SOLID) {
            double value = - factor * (vf.get_v(i, j) - bc.prescribed_v(i, j));
            poisson_rhs(cell_id) += value;
        }

        if (bc.type(i, j + 1) == CellType::FLUID) {
            poisson_matrix(cell_id, cell_id) ++;
            poisson_matrix(cell_id, convert_index(i, j + 1)) --;
        }
        else if (bc.type(i, j + 1) == CellType::BOUNDARY) {
            poisson_matrix(cell_id, cell_id) ++;
            poisson_rhs(cell_id) += pf.outside_pressure();
        }
        else if (bc.type(i, j + 1) == CellType::SOLID) {
            double value = factor * (vf.get_v(i, j + 1) - bc.prescribed_v(i, j + 1));
            poisson_rhs(cell_id) += value;
        }
    }

    void Simulator::apply_pressure_gradient(double timestep) {
        std::size_t w = grid_.width(), h = grid_.height();
        const auto& pf = grid_.pressure();
        auto& vf = grid_.velocity();

        const double factor = timestep / (fluid_density_ * grid_.resolution());

        /// Apply pressure gradient across vertical faces
        for (int i = 0; i <= w; i ++) 
            for (int j = 0; j < h; j ++) {
                double delta_pressure = 
                    pf.read_p_or_outside(i, j) - pf.read_p_or_outside(i - 1, j);
                vf.get_u(i, j) -= factor * delta_pressure;
            }
        
        /// Apply pressure gradient across horizontal faces
        for (int i = 0; i < w; i ++) 
            for (int j = 0; j <= h; j ++) {
                double delta_pressure = 
                    pf.read_p_or_outside(i, j) - pf.read_p_or_outside(i, j - 1);
                vf.get_v(i, j) -= factor * delta_pressure;
            }
    }

    void Simulator::set_pressure_values(const linalg::Vector& pressure_values) {
        int w = grid_.width(), h = grid_.height();
        auto& pf = grid_.pressure();

        for (int i = 0; i < w; i ++)
            for (int j = 0; j < h; j ++) {
                std::size_t id = static_cast<std::size_t>(i) * h + 
                                    static_cast<std::size_t>(j);
                pf.get_p(i, j) = pressure_values(id);
            }
    }

    void Simulator::project(double timestep) {
        /* Determine the pressure in all cells of the grid to 
         * make the velocity field divergence free, taking
         * into account the boundary condtions.
         */

        int number_of_cells = grid_.width() * grid_.height();
        int w = grid_.width(), h = grid_.height();
        const auto& bc = grid_.boundary_conditions();

        linalg::Matrix poisson_matrix(           // Each connected island of fluid cells requires
            number_of_cells + 1, number_of_cells // an additional equation in the system (Neumann problem)
        );                                       // Currently, only one island is considered

        linalg::Vector poisson_rhs(number_of_cells + 1);
        
        auto is_edge = [&](int i, int j) { 
            return (bc.type(i - 1, j) == CellType::BOUNDARY ||
                    bc.type(i, j - 1) == CellType::BOUNDARY ||
                    bc.type(i + 1, j) == CellType::BOUNDARY ||
                    bc.type(i, j + 1) == CellType::BOUNDARY);
        };

        bool is_island = true;
        for (int i = 0; i < w; i ++)
            for (int j = 0; j < h; j ++) {
                build_equation(i, j, timestep, poisson_matrix, poisson_rhs);

                if (bc.type(i, j) == CellType::FLUID && is_edge(i, j))
                    is_island = false;
            }
        
        if (is_island) {
            /// We add an additional equation (mean pressure value to be 0) to fix solution
            poisson_rhs(number_of_cells) = 0;
            for (int i = 0; i < number_of_cells; i ++) // Should apply only for fluid cells but 
                poisson_matrix(number_of_cells, i) = 1;        // Solid cells have zero pressure
        }

        linalg::Vector pressure_values = linalg::gaussian_elimination(
            poisson_matrix, poisson_rhs
        );

        set_pressure_values(pressure_values);
        apply_pressure_gradient(timestep);
    }

    void Simulator::tick() {
        double timestep = determine_timestep();
        
        advect(timestep);
        apply_boundary_conditions();

        if (apply_gravity_) {
            apply_body_forces(timestep, 0.0, GRAVITY_AY);
            apply_boundary_conditions();
        }

        project(timestep);
        //apply_boundary_conditions();
    }

    Grid& Simulator::grid() {
        return grid_;
    }

    const Grid& Simulator::grid() const {
        return grid_;
    }
}