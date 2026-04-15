#include "core/velocity_field.hpp"
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <algorithm>

namespace cfd {

    VelocityField::VelocityField(
        std::size_t width, 
        std::size_t height, 
        double resolution
    )
        : width_(width), 
          height_(height), 
          resolution_(resolution), 
          u_((width + 1) * height), 
          v_(width * (height + 1)),
          outside_velocity_({0, 0}) {
    }
    
    
    void VelocityField::check_u_bounds(int i, int j) const {
        /// Check if the indexes are inside the u field bounds

        const int width = static_cast<int>(width_);
        const int height = static_cast<int>(height_);

        if (i < 0 || i > width || j < 0 || j >= height) {
            throw std::out_of_range("u field index out of bounds");
        }
    }

    void VelocityField::check_v_bounds(int i, int j) const {
        /// Check if the indexes are inside the v field bounds

        const int width = static_cast<int>(width_);
        const int height = static_cast<int>(height_);

        if (i < 0 || i >= width || j < 0 || j > height) {
            throw std::out_of_range("v field index out of bounds");
        }
    }

    void VelocityField::check_coordinates_bounds(double x, double y) const {
        const double max_x = static_cast<double>(width_) * resolution_;
        const double max_y = static_cast<double>(height_) * resolution_;

        if (x < -EPS || x > max_x + EPS ||
            y < -EPS || y > max_y + EPS) {
            throw std::out_of_range("coordinates out of bounds");
        }
    }

    double& VelocityField::get_u(int i, int j) {
        // Return the u-velocity at the left face of cell [i, j].

        check_u_bounds(i, j);
        return u_[i * height_ + j];
    }

    const double& VelocityField::get_u(int i, int j) const {
        // Return the u-velocity at the left face of cell [i, j].

        check_u_bounds(i, j);
        return u_[i * height_ + j];
    }

    double& VelocityField::get_v(int i, int j) {
        // Return the v-velocity at the bottom face of cell [i, j].

        check_v_bounds(i, j);
        return v_[i * (height_ + 1) + j];
    }

    const double& VelocityField::get_v(int i, int j) const {
        // Return the v-velocity at the bottom face of cell [i, j].

        check_v_bounds(i, j);
        return v_[i * (height_ + 1) + j];
    }


    double VelocityField::read_u_or_outside(int i, int j) const {
        /// Returns the u component in cell [i, j]

        const int width = static_cast<int>(width_);
        const int height = static_cast<int>(height_);
        
        /// Clamp to the closest cell if indexes are outside bounds
        
        i = std::clamp(i, 0, width);
        j = std::clamp(j, 0, height - 1);

        return u_[i * height_ + j];
    }

    double VelocityField::read_v_or_outside(int i, int j) const {
        /// Returns the v component in cell [i, j]

        const int width = static_cast<int>(width_);
        const int height = static_cast<int>(height_);

        /// Clamp to the closest cell if indexes are outside bounds

        i = std::clamp(i, 0, width - 1); 
        j = std::clamp(j, 0, height);

        return v_[i * (height_ + 1) + j];
    }

    double VelocityField::get_max_velocity_component() const {
        /// Return the maximum velocity component in the velocity field

        double max_vel = 0.0;

        // Vertical faces.
        for (std::size_t i = 0; i <= width_; ++i) {
            for (std::size_t j = 0; j < height_; ++j) {
                const double u = std::abs(get_u(i, j));
                if (u > max_vel) {
                    max_vel = u;
                }
            }
        }

        // Horizontal faces.
        for (std::size_t i = 0; i < width_; ++i) {
            for (std::size_t j = 0; j <= height_; ++j) {
                const double v = std::abs(get_v(i, j));
                if (v > max_vel) {
                    max_vel = v;
                }
            }
        }

        return max_vel;
    }

    void VelocityField::set_outside_velocity(Velocity2D outside_velocity) {
        outside_velocity_ = outside_velocity;
    }
    
    double VelocityField::get_divergence(int i, int j) const {
        // Return the Divergence of the cell [i, j]

        const double du_dx = (get_u(i + 1, j) - get_u(i, j)) / resolution_;
        const double dv_dy = (get_v(i, j + 1) - get_v(i, j)) / resolution_;

        return du_dx + dv_dy;
    }

    double VelocityField::get_u_laplacian(int i, int j) const {
        // Return the Laplacian of u at the left face of cell [i, j].

        const double u_center = read_u_or_outside(i, j);

        const double x_contribution =
            read_u_or_outside(i + 1, j) +
            read_u_or_outside(i - 1, j) -
            2.0 * u_center;

        const double y_contribution =
            read_u_or_outside(i, j + 1) +
            read_u_or_outside(i, j - 1) -
            2.0 * u_center;

        return (x_contribution + y_contribution) / (resolution_ * resolution_);
    }

    double VelocityField::get_v_laplacian(int i, int j) const {
        // Return the Laplacian of v at the bottom face of cell [i, j].

        const double v_center = read_v_or_outside(i, j);

        const double x_contribution =
            read_v_or_outside(i + 1, j) +
            read_v_or_outside(i - 1, j) -
            2.0 * v_center;

        const double y_contribution =
            read_v_or_outside(i, j + 1) +
            read_v_or_outside(i, j - 1) -
            2.0 * v_center;

        return (x_contribution + y_contribution) / (resolution_ * resolution_);
    }

    Velocity2D VelocityField::sample_at_vertical_face(int i, int j) const {
        // Return the velocity at the vertical face on the left of cell [i, j].

        const double u = read_u_or_outside(i, j);

        const double v =
            (read_v_or_outside(i - 1, j) +
            read_v_or_outside(i - 1, j + 1) +
            read_v_or_outside(i, j) +
            read_v_or_outside(i, j + 1)) /
            4.0;

        return {u, v};
    }

    Velocity2D VelocityField::sample_at_horizontal_face(int i, int j) const {
        // Return the velocity at the horizontal face on the bottom of cell [i, j].

        const double u =
            (read_u_or_outside(i, j) +
            read_u_or_outside(i + 1, j) +
            read_u_or_outside(i, j - 1) +
            read_u_or_outside(i + 1, j - 1)) /
            4.0;

        const double v = read_v_or_outside(i, j);

        return {u, v};
    }

    Velocity2D VelocityField::sample_at_corner(int i, int j) const {
        // Return the velocity at the bottom-left corner of cell [i, j].

        const Velocity2D up = sample_at_vertical_face(i, j);
        const Velocity2D right = sample_at_horizontal_face(i, j);
        const Velocity2D down = sample_at_vertical_face(i, j - 1);
        const Velocity2D left = sample_at_horizontal_face(i - 1, j);

        const double u = (up.u + down.u) / 2.0;
        const double v = (left.v + right.v) / 2.0;

        return {u, v};
    }

    Velocity2D VelocityField::sample_at_coordinates(double x, double y) const {
        // Return the velocity at the given coordinates using bilinear interpolation.

        const int i = static_cast<int>(std::floor(x / resolution_));
        const int j = static_cast<int>(std::floor(y / resolution_));

        const double alpha_x =
            (x - static_cast<double>(i) * resolution_) / resolution_;
        const double alpha_y =
            (y - static_cast<double>(j) * resolution_) / resolution_;

        const Velocity2D ld = sample_at_corner(i, j);
        const Velocity2D lu = sample_at_corner(i, j + 1);
        const Velocity2D ru = sample_at_corner(i + 1, j + 1);
        const Velocity2D rd = sample_at_corner(i + 1, j);

        const double u =
            (1.0 - alpha_x) * (1.0 - alpha_y) * ld.u +
            alpha_x * (1.0 - alpha_y) * rd.u +
            (1.0 - alpha_x) * alpha_y * lu.u +
            alpha_x * alpha_y * ru.u;

        const double v =
            (1.0 - alpha_x) * (1.0 - alpha_y) * ld.v +
            alpha_x * (1.0 - alpha_y) * rd.v +
            (1.0 - alpha_x) * alpha_y * lu.v +
            alpha_x * alpha_y * ru.v;

        return {u, v};
    }

}