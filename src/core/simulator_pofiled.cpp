#include "core/simulator.hpp"
#include "linalg/gaussian_elimination.hpp"
#include "linalg/conjugate_gradient.hpp"
#include <optional>
#include <stdexcept>
#include <iostream>
#include <chrono>
#include <iomanip>

namespace cfd {

    Simulator::Simulator(std::size_t width,
                     std::size_t height,
                     double resolution,
                     double fluid_density,
                     bool apply_gravity,
                     double viscosity)
    : grid_(width, height, resolution),
      fluid_density_(fluid_density),
      apply_gravity_(apply_gravity),
      viscosity_(viscosity) {}

    double Simulator::determine_timestep() const {
        /// Determine timestep using CFL condition

        const double max_velocity = grid_.velocity().get_max_velocity_component();
        const double dx = grid_.resolution();

        double advective_dt = FALLBACK_MAX_TIMESTEP;
        if (max_velocity > 0.0) {
            advective_dt = CFL * dx / max_velocity;
        }

        double diffusive_dt = FALLBACK_MAX_TIMESTEP;
        if (viscosity_ > 0.0) {
            diffusive_dt = 0.25 * dx * dx / viscosity_;
        }

        return std::min(
            FALLBACK_MAX_TIMESTEP,
            std::min(advective_dt, diffusive_dt)
        );
    }

    void Simulator::advect(double timestep) {
        // Perform semi-Lagrangian advection on the velocity field.

        const auto& old_velocity = grid_.velocity();
        const std::size_t width = grid_.width();
        const std::size_t height = grid_.height();
        const double dx = grid_.resolution();

        VelocityField next_velocity(width, height, dx);

        // Advect vertical faces.
        for (std::size_t i = 0; i <= width; i++) {
            for (std::size_t j = 0; j < height; j++) {
                const auto velocity = old_velocity.sample_at_vertical_face(i, j);

                const double x = static_cast<double>(i) * dx;
                const double y = (static_cast<double>(j) + 0.5) * dx;

                const double previous_x = x - timestep * velocity.u;
                const double previous_y = y - timestep * velocity.v;

                const auto advected_velocity =
                    old_velocity.sample_at_coordinates(previous_x, previous_y);

                next_velocity.get_u(i, j) = advected_velocity.u;
            }
        }

        // Advect horizontal faces.
        for (std::size_t i = 0; i < width; i++) {
            for (std::size_t j = 0; j <= height; j++) {
                const auto velocity = old_velocity.sample_at_horizontal_face(i, j);

                const double x = (static_cast<double>(i) + 0.5) * dx;
                const double y = static_cast<double>(j) * dx;

                const double previous_x = x - timestep * velocity.u;
                const double previous_y = y - timestep * velocity.v;

                const auto advected_velocity =
                    old_velocity.sample_at_coordinates(previous_x, previous_y);

                next_velocity.get_v(i, j) = advected_velocity.v;
            }
        }

        grid_.velocity() = std::move(next_velocity);
    }

    void Simulator::apply_body_forces(double timestep, double ax, double ay) {
        /// Apply body forces on the velocity field

        auto& velocity = grid_.velocity();
        const std::size_t width = grid_.width();
        const std::size_t height = grid_.height();

        // Apply on vertical faces.
        for (std::size_t i = 0; i <= width; i++) {
            for (std::size_t j = 0; j < height; j++) {
                velocity.get_u(i, j) += timestep * ax;
            }
        }

        // Apply on horizontal faces.
        for (std::size_t i = 0; i < width; i++) {
            for (std::size_t j = 0; j <= height; j++) {
                velocity.get_v(i, j) += timestep * ay;
            }
        }
    }

    void Simulator::apply_boundary_conditions() {
        /// Apply boundary conditions on the velocity field

        const auto& boundary_conditions = grid_.boundary_conditions();
        auto& velocity = grid_.velocity();

        const auto& prescribed_u = boundary_conditions.prescribed_u();
        const auto& prescribed_v = boundary_conditions.prescribed_v();

        const std::size_t height = grid_.height();

        for (const auto& [key, value] : prescribed_u) {
            const std::size_t i = key / height;
            const std::size_t j = key % height;

            velocity.get_u(i, j) = value;
        }

        for (const auto& [key, value] : prescribed_v) {
            const std::size_t i = key / (height + 1);
            const std::size_t j = key % (height + 1);

            velocity.get_v(i, j) = value;
        }
    }

    std::optional<double> Simulator::get_boundary_pressure_constraint(
        int i,
        int j,
        double timestep
    ) {
        // If a face between a boundary and a cell has a prescribed velocity
        // it imposes a fixed pressure constraint on the cell
        
        const auto& bc = grid_.boundary_conditions();
        const auto& vf = grid_.velocity();

        std::optional<double> constraint;

        auto add_pressure_constraint = [&](double pressure) {
            if (!constraint.has_value()) {
                constraint = pressure;
                return;
            }

            if (std::abs(*constraint - pressure) > TOL) {
                throw std::runtime_error(
                    "Incompatible pressure/boundary conditions:\n"
                    "cell pressure is overconstrained"
                );
            }
        };

        const double factor = fluid_density_ * grid_.resolution() / timestep;

        if (bc.type(i - 1, j) == CellType::BOUNDARY &&
            bc.is_u_prescribed(i, j)) {
            const double cell_pressure =
                bc.prescribed_p(i - 1, j) +
                factor * (vf.get_u(i, j) - bc.prescribed_u(i, j));

            add_pressure_constraint(cell_pressure);
        }

        if (bc.type(i, j - 1) == CellType::BOUNDARY &&
            bc.is_v_prescribed(i, j)) {
            const double cell_pressure =
                bc.prescribed_p(i, j - 1) +
                factor * (vf.get_v(i, j) - bc.prescribed_v(i, j));

            add_pressure_constraint(cell_pressure);
        }

        if (bc.type(i + 1, j) == CellType::BOUNDARY &&
            bc.is_u_prescribed(i + 1, j)) {
            const double cell_pressure =
                bc.prescribed_p(i + 1, j) -
                factor * (vf.get_u(i + 1, j) - bc.prescribed_u(i + 1, j));

            add_pressure_constraint(cell_pressure);
        }

        if (bc.type(i, j + 1) == CellType::BOUNDARY &&
            bc.is_v_prescribed(i, j + 1)) {
            const double cell_pressure =
                bc.prescribed_p(i, j + 1) -
                factor * (vf.get_v(i, j + 1) - bc.prescribed_v(i, j + 1));

            add_pressure_constraint(cell_pressure);
        }

        return constraint;
    }

    double Simulator::get_pressure(int i, int j) {
        // Return the pressure in cell [i, j]
        // The cell could be of type boundary

        int width = static_cast<int>(grid_.width());
        int height = static_cast<int>(grid_.height());

        const auto& bc = grid_.boundary_conditions();
        const auto& pf = grid_.pressure();

        if (i < 0 || i >= width ||
            j < 0 || j >= height) {
                return bc.prescribed_p(i, j); 
            }
        
        return pf.get_p(i, j);
    }
    

    void fix_cell(
        int cell_id,
        double val,
        int num_of_cells,
        linalg::LinearOperator& poisson_matrix,
        linalg::Vector& poisson_rhs
    ) {
        // Fixes the pressure of a cell to a given value
        // in the poisson set of equations
        
        for (int row = 0; row < num_of_cells; row++) {
            if (row != cell_id) {
                poisson_rhs(row) -= poisson_matrix(row, cell_id) * val;
            }
        }

        for (int id = 0; id < num_of_cells; id++) {
            poisson_matrix(cell_id, id) = 0;
            poisson_matrix(id, cell_id) = 0;
        }

        poisson_matrix(cell_id, cell_id) = 1;
        poisson_rhs(cell_id) = val;
    }

    std::optional<double> Simulator::build_equation(
        int i,
        int j,
        double timestep,
        linalg::LinearOperator& poisson_matrix,
        linalg::Vector& poisson_rhs
    ) {
        // Build the equation for cell [i, j]
        
        const std::size_t h = grid_.height();

        auto cell_index = [h](std::size_t i, std::size_t j) {
            return i * h + j;
        };

        const auto& bc = grid_.boundary_conditions();
        const auto& vf = grid_.velocity();

        const std::size_t cell_id =
            cell_index(static_cast<std::size_t>(i), static_cast<std::size_t>(j));

        if (bc.type(i, j) == CellType::SOLID) {
            // Enforce zero pressure on solid cells.
            poisson_matrix(cell_id, cell_id) = 1;
            poisson_rhs(cell_id) = 0;
            return std::nullopt;
        }

        if (bc.is_p_prescribed(i, j)) {
            return bc.prescribed_p(i, j);
        }

        const std::optional<double> boundary_constraint =
            get_boundary_pressure_constraint(i, j, timestep);

        if (boundary_constraint.has_value()) {
            return *boundary_constraint;
        }

        const double div = vf.get_divergence(i, j);

        poisson_rhs(cell_id) =
            -(fluid_density_ * grid_.resolution() * grid_.resolution() * div) /
            timestep;

        const double factor = fluid_density_ * grid_.resolution() / timestep;

        if (bc.type(i - 1, j) == CellType::FLUID) {
            poisson_matrix(cell_id, cell_id)++;
            poisson_matrix(cell_id, cell_index(i - 1, j))--;
        } else if (bc.type(i - 1, j) == CellType::BOUNDARY) {
            poisson_matrix(cell_id, cell_id)++;
            poisson_rhs(cell_id) += bc.prescribed_p(i - 1, j);
        } else if (bc.type(i - 1, j) == CellType::SOLID) {
            const double value = -factor * (vf.get_u(i, j) - bc.prescribed_u(i, j));
            poisson_rhs(cell_id) += value;
        }

        if (bc.type(i + 1, j) == CellType::FLUID) {
            poisson_matrix(cell_id, cell_id)++;
            poisson_matrix(cell_id, cell_index(i + 1, j))--;
        } else if (bc.type(i + 1, j) == CellType::BOUNDARY) {
            poisson_matrix(cell_id, cell_id)++;
            poisson_rhs(cell_id) += bc.prescribed_p(i + 1, j);
        } else if (bc.type(i + 1, j) == CellType::SOLID) {
            const double value = factor * (vf.get_u(i + 1, j) - bc.prescribed_u(i + 1, j));
            poisson_rhs(cell_id) += value;
        }

        if (bc.type(i, j - 1) == CellType::FLUID) {
            poisson_matrix(cell_id, cell_id)++;
            poisson_matrix(cell_id, cell_index(i, j - 1))--;
        } else if (bc.type(i, j - 1) == CellType::BOUNDARY) {
            poisson_matrix(cell_id, cell_id)++;
            poisson_rhs(cell_id) += bc.prescribed_p(i, j - 1);
        } else if (bc.type(i, j - 1) == CellType::SOLID) {
            const double value = -factor * (vf.get_v(i, j) - bc.prescribed_v(i, j));
            poisson_rhs(cell_id) += value;
        }

        if (bc.type(i, j + 1) == CellType::FLUID) {
            poisson_matrix(cell_id, cell_id)++;
            poisson_matrix(cell_id, cell_index(i, j + 1))--;
        } else if (bc.type(i, j + 1) == CellType::BOUNDARY) {
            poisson_matrix(cell_id, cell_id)++;
            poisson_rhs(cell_id) += bc.prescribed_p(i, j + 1);
        } else if (bc.type(i, j + 1) == CellType::SOLID) {
            const double value = factor * (vf.get_v(i, j + 1) - bc.prescribed_v(i, j + 1));
            poisson_rhs(cell_id) += value;
        }

        return std::nullopt;
    }

    void Simulator::apply_pressure_gradient(double timestep) {
        const std::size_t w = grid_.width();
        const std::size_t h = grid_.height();
        auto& vf = grid_.velocity();

        const double factor = timestep / (fluid_density_ * grid_.resolution());

        // Apply pressure gradient across vertical faces.
        for (int i = 0; i <= w; i++) {
            for (int j = 0; j < h; j++) {
                const double delta_pressure =
                    get_pressure(i, j) - get_pressure(i - 1, j);

                vf.get_u(i, j) -= factor * delta_pressure;
            }
        }

        // Apply pressure gradient across horizontal faces.
        for (int i = 0; i < w; i++) {
            for (int j = 0; j <= h; j++) {
                const double delta_pressure =
                    get_pressure(i, j) - get_pressure(i, j - 1);

                vf.get_v(i, j) -= factor * delta_pressure;
            }
        }
    }

    void Simulator::set_pressure_values(const linalg::Vector& pressure_values) {
        const int w = grid_.width();
        const int h = grid_.height();
        auto& pf = grid_.pressure();

        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {
                const std::size_t id =
                    static_cast<std::size_t>(i) * h +
                    static_cast<std::size_t>(j);

                pf.get_p(i, j) = pressure_values(id);
            }
        }
    }

    void Simulator::project(double timestep) {
        using clock = std::chrono::steady_clock;

        auto elapsed_ms = [](const clock::time_point& a, const clock::time_point& b) {
            return std::chrono::duration<double, std::milli>(b - a).count();
        };

        auto measure = [&](double& total_ms, int& calls, auto&& fn) {
            auto t0 = clock::now();
            fn();
            auto t1 = clock::now();
            total_ms += elapsed_ms(t0, t1);
            ++calls;
        };

        double setup_ms = 0.0;
        double main_loop_ms = 0.0;
        double build_equation_ms = 0.0;
        double island_fix_ms = 0.0;
        double constrained_fix_ms = 0.0;
        double sterilize_ms = 0.0;
        double cg_ms = 0.0;
        double set_pressure_ms = 0.0;
        double apply_gradient_ms = 0.0;

        int setup_calls = 0;
        int main_loop_calls = 0;
        int build_equation_calls = 0;
        int island_fix_calls = 0;
        int constrained_fix_calls = 0;
        int sterilize_calls = 0;
        int cg_calls = 0;
        int set_pressure_calls = 0;
        int apply_gradient_calls = 0;

        auto project_start = clock::now();

        int number_of_cells = 0;
        int w = 0;
        int h = 0;

        const BoundaryConditions* bc_ptr = nullptr;

        linalg::LinearOperator poisson_matrix(0, 0);
        linalg::Vector poisson_rhs(0);
        std::vector<std::pair<int, double>> constrained_cells;

        measure(setup_ms, setup_calls, [&] {
            number_of_cells = grid_.width() * grid_.height();
            w = grid_.width();
            h = grid_.height();
            bc_ptr = &grid_.boundary_conditions();

            poisson_matrix = linalg::LinearOperator(number_of_cells, number_of_cells);
            poisson_rhs = linalg::Vector(number_of_cells);
            constrained_cells.clear();
        });

        const auto& bc = *bc_ptr;

        auto is_edge = [&](int i, int j) {
            return (
                bc.type(i - 1, j) == CellType::BOUNDARY ||
                bc.type(i, j - 1) == CellType::BOUNDARY ||
                bc.type(i + 1, j) == CellType::BOUNDARY ||
                bc.type(i, j + 1) == CellType::BOUNDARY
            );
        };

        bool is_island = true;

        measure(main_loop_ms, main_loop_calls, [&] {
            for (int i = 0; i < w; i++) {
                for (int j = 0; j < h; j++) {
                    std::optional<double> constraint;

                    measure(build_equation_ms, build_equation_calls, [&] {
                        constraint = build_equation(i, j, timestep, poisson_matrix, poisson_rhs);
                    });

                    if (bc.type(i, j) == CellType::FLUID && is_edge(i, j)) {
                        is_island = false;
                    }

                    if (constraint.has_value()) {
                        is_island = false;
                        const int cell_id = i * h + j;
                        constrained_cells.push_back({cell_id, *constraint});
                    }
                }
            }
        });

        if (is_island) {
            measure(island_fix_ms, island_fix_calls, [&] {
                bool found = false;

                for (int i = 0; i < w && !found; i++) {
                    for (int j = 0; j < h && !found; j++) {
                        if (bc.type(i, j) == CellType::FLUID) {
                            const int cell_id = i * h + j;
                            constrained_cells.push_back({cell_id, 0.0});

                            fix_cell(cell_id, 0.0, number_of_cells, poisson_matrix, poisson_rhs);
                            found = true;
                        }
                    }
                }
            });
        }

        measure(constrained_fix_ms, constrained_fix_calls, [&] {
            for (const auto& [id, value] : constrained_cells) {
                fix_cell(id, value, number_of_cells, poisson_matrix, poisson_rhs);
            }
        });

        measure(sterilize_ms, sterilize_calls, [&] {
            poisson_matrix.sterilize();
        });

        linalg::Vector pressure_values(number_of_cells);

        measure(cg_ms, cg_calls, [&] {
            pressure_values = linalg::conjugate_gradient(poisson_matrix, poisson_rhs);
        });

        measure(set_pressure_ms, set_pressure_calls, [&] {
            set_pressure_values(pressure_values);
        });

        measure(apply_gradient_ms, apply_gradient_calls, [&] {
            apply_pressure_gradient(timestep);
        });

        auto project_end = clock::now();
        double total_project_ms = elapsed_ms(project_start, project_end);

        auto print_stat = [&](const char* name, double ms, int calls) {
            double pct = (total_project_ms > 0.0) ? (100.0 * ms / total_project_ms) : 0.0;

            std::cout
                << std::left << std::setw(28) << name
                << " | total: " << std::setw(10) << ms << " ms"
                << " | pct: " << std::setw(8) << pct << " %"
                << " | calls: " << std::setw(6) << calls;

            if (calls > 1) {
                std::cout << " | avg/call: " << (ms / calls) << " ms";
            }

            std::cout << '\n';
        };

        std::cout << std::fixed << std::setprecision(3);
        std::cout << "\n=== Project Profile ===\n";
        std::cout << "Total project time: " << total_project_ms << " ms\n";

        print_stat("setup", setup_ms, setup_calls);
        print_stat("main grid loop", main_loop_ms, main_loop_calls);
        print_stat("build_equation", build_equation_ms, build_equation_calls);

        if (island_fix_calls > 0) {
            print_stat("island fix", island_fix_ms, island_fix_calls);
        }

        print_stat("fix constrained cells", constrained_fix_ms, constrained_fix_calls);
        print_stat("sterilize", sterilize_ms, sterilize_calls);
        print_stat("conjugate_gradient", cg_ms, cg_calls);
        print_stat("set_pressure_values", set_pressure_ms, set_pressure_calls);
        print_stat("apply_pressure_gradient", apply_gradient_ms, apply_gradient_calls);
        std::cout << "=======================\n";

        set_pressure_values(pressure_values);
        apply_pressure_gradient(timestep);
    }

    void Simulator::tick() {
        
        const double timestep = determine_timestep();
        
        advect(timestep);
        apply_boundary_conditions();
        
        apply_viscosity(timestep);
        apply_boundary_conditions();
        
        if (apply_gravity_) {
            apply_body_forces(timestep, 0.0, GRAVITY_AY);
            apply_boundary_conditions();
        }
        
        project(timestep);
        apply_boundary_conditions();
    }

    void Simulator::apply_viscosity(double timestep) {
        const std::size_t h = grid_.height();
        const std::size_t w = grid_.width();

        auto& vf = grid_.velocity();
        const double dx = grid_.resolution();

        VelocityField next_field(w, h, dx);

        // Apply viscosity on vertical faces.
        for (int i = 0; i <= w; i++)
            for (int j = 0; j < h; j++) {
                next_field.get_u(i, j) =
                    vf.get_u(i, j) + timestep * viscosity_ * vf.get_u_laplacian(i, j);
            }

        // Apply viscosity on horizontal faces.
        for (int i = 0; i < w; i++)
            for (int j = 0; j <= h; j++) {
                next_field.get_v(i, j) =
                    vf.get_v(i, j) + timestep * viscosity_ * vf.get_v_laplacian(i, j);
            }

        vf = std::move(next_field);
    }

    Grid& Simulator::grid() {
        return grid_;
    }

    const Grid& Simulator::grid() const {
        return grid_;
    }
}