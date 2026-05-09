#include "setups/airfoil_flow.hpp"

#include <cmath>
#include <stdexcept>

namespace cfd::setups {

    AirfoilFlowSetup::AirfoilFlowSetup(
        const std::string& naca_code,
        double chord,
        double freestream_velocity,
        double dx,
        double angle_of_attack_deg
    )
        : airfoil_(naca_code, chord),
          v_air_(freestream_velocity),
          dx_(dx),
          aoa_deg_(angle_of_attack_deg)
    {
        if (dx_ <= 0.0) {
            throw std::invalid_argument("dx must be positive");
        }

        if (v_air_ < 0.0) {
            throw std::invalid_argument("freestream velocity cannot be negative");
        }
    }

    void AirfoilFlowSetup::impose_boundary_conditions(Grid& grid) const
    {
        auto& bc = grid.boundary_conditions();
        auto& vf = grid.velocity();

        const int w = grid.width();
        const int h = grid.height();

        // --------------------------------------------------
        // Initialize whole velocity field to freestream
        // --------------------------------------------------
        for (int j = 0; j < h; ++j) {
            for (int i = 0; i <= w; ++i) {
                vf.set_u(i, j, v_air_);
            }
        }

        for (int j = 0; j <= h; ++j) {
            for (int i = 0; i < w; ++i) {
                vf.set_v(i, j, 0.0);
            }
        }

        // --------------------------------------------------
        // Left inlet velocity
        // --------------------------------------------------
        for (int j = 0; j < h; ++j) {
            bc.prescribe_u_value(0, j, v_air_);
        }

        // --------------------------------------------------
        // Top and bottom far-field
        // --------------------------------------------------
        for (int i = 0; i < w; ++i) {
            bc.prescribe_v_value(i, 0, 0.0);
            bc.prescribe_v_value(i, h, 0.0);

            bc.prescribe_u_value(i, 0, v_air_);
            bc.prescribe_u_value(i, h - 1, v_air_);
        }

        // --------------------------------------------------
        // Right pressure outlet
        // --------------------------------------------------
        for (int j = 0; j < h; ++j) {
            bc.prescribe_p_value(w, j, 0.0);
            bc.prescribe_p_value(-1, j, 0.0);
        }

        // --------------------------------------------------
        // Airfoil solid cells
        // --------------------------------------------------
        for (int j = 0; j < h; ++j) {
            for (int i = 0; i < w; ++i) {
                if (is_inside_airfoil(i, j, w, h)) {
                    make_solid_cell(i, j, grid, bc);
                }
            }
        }
    }

    bool AirfoilFlowSetup::is_inside_airfoil(
        int i,
        int j,
        int w,
        int h
    ) const
    {
        constexpr double pi = 3.14159265358979323846;

        const double domain_width = w * dx_;
        const double domain_height = h * dx_;

        const double alpha = -aoa_deg_ * pi / 180.0;

        const double x_le = 0.5 * domain_width - 0.5 * airfoil_.chord();
        const double y_chord = 0.5 * domain_height;

        // Rotate around quarter-chord point
        const double x_pivot = x_le + 0.25 * airfoil_.chord();
        const double y_pivot = y_chord;

        // Cell center in world coordinates
        const double x = (static_cast<double>(i) + 0.5) * dx_;
        const double y = (static_cast<double>(j) + 0.5) * dx_;

        // Translate relative to pivot
        const double xr = x - x_pivot;
        const double yr = y - y_pivot;

        // Inverse-rotate world point into airfoil-local coordinates
        const double x_local =
            x_pivot
            + xr * std::cos(alpha)
            + yr * std::sin(alpha);

        const double y_local =
            y_pivot
            - xr * std::sin(alpha)
            + yr * std::cos(alpha);

        const double X = (x_local - x_le) / airfoil_.chord();

        if (X < 0.0 || X > 1.0) {
            return false;
        }

        const double y_c = airfoil_.mean_camber_line(X);
        const double y_t = airfoil_.half_thickness(X);

        const double y_airfoil_center = y_chord + y_c;

        return (y_airfoil_center - y_t <= y_local)
            && (y_local <= y_airfoil_center + y_t);
    }

    void AirfoilFlowSetup::make_solid_cell(
        int i,
        int j,
        Grid& grid,
        BoundaryConditions& bc
    ) const
    {
        auto& vf = grid.velocity();

        set_to_solid(i, j, bc);

        vf.set_u(i, j, 0.0);
        vf.set_u(i + 1, j, 0.0);

        vf.set_v(i, j, 0.0);
        vf.set_v(i, j + 1, 0.0);

        bc.prescribe_u_value(i, j, 0.0);
        bc.prescribe_u_value(i + 1, j, 0.0);

        bc.prescribe_v_value(i, j, 0.0);
        bc.prescribe_v_value(i, j + 1, 0.0);
    }

}
