#include "core/velocity_field.hpp"
#include <stdexcept>
#include <cmath>
#include <iostream>

namespace cfd {

    VelocityField::VelocityField(
        std::size_t width, std::size_t height, double resolution
    ): width_(width), height_(height), resolution_(resolution), 
        u_((width + 1) * height), v_(width * (height + 1)),
        outside_velocity_({0, 0}) {}
    
    
    void VelocityField::check_u_bounds(int i, int j) const {
        if (i > static_cast<int>(width_) || j >= static_cast<int>(height_) 
            || i < 0 || j < 0) {
            throw std::out_of_range("u field index out of bounds");
        }
    }
    
    void VelocityField::check_v_bounds(int i, int j) const {
        if (i >= static_cast<int>(width_) || j > static_cast<int>(height_) 
            || i < 0 || j < 0) {
            throw std::out_of_range("v field index out of bounds");
        }
    }

    void VelocityField::check_coordinates_bounds(double x, double y) const {
        if (x < -EPS || x > static_cast<double>(width_) * resolution_ + EPS) {
            throw std::out_of_range("coordinates out of bounds");
        }
     
        if (y < -EPS || y > static_cast<double>(height_) * resolution_ + EPS) {
            throw std::out_of_range("coordinates out of bounds");
        }
    }

    double& VelocityField::get_u(int i, int j) {
        /// get u component of velocity at the left face of cell [i, j]
        check_u_bounds(i, j);
        return u_[i * height_ + j];
    }

    const double& VelocityField::get_u(int i, int j) const {
        /// get u component of velocity at the left face of cell [i, j]
        check_u_bounds(i, j);
        return u_[i * height_ + j];
    }

    double VelocityField::read_u_or_outside(int i, int j) const {
        if (i > static_cast<int>(width_) || j >= static_cast<int>(height_) 
            || i < 0 || j < 0)
            return outside_velocity_.u;

        return u_[i * height_ + j];
    }

    double VelocityField::read_v_or_outside(int i, int j) const {
        if (i >= static_cast<int>(width_) || j > static_cast<int>(height_) 
            || i < 0 || j < 0)
            return outside_velocity_.v;
        
        return v_[i * (height_ + 1) + j];
    }

    double& VelocityField::get_v(int i, int j) {
        /// get v component of velocity at the bottom face of cell [i, j]
        check_v_bounds(i, j);
        return v_[i * (height_ + 1) + j];
    }

    const double& VelocityField::get_v(int i, int j) const {
        /// get v component of velocity at the bottom face of cell [i, j]
        check_v_bounds(i, j);
        return v_[i * (height_ + 1) + j];
    }

    Velocity2D VelocityField::sample_at_vertical_face(int i, int j) const {
        /// Get velocity vector at the vertical face on the left of cell [i, j]
        
        /// interpolates either the 2 or 4 closes horizontal faces to obtain v
        double v = (read_v_or_outside(i - 1, j) + read_v_or_outside(i - 1, j + 1) +
            read_v_or_outside(i, j) + read_v_or_outside(i, j + 1)) / 4.0;
        double u = read_u_or_outside(i ,j);

        return {u, v};
    }

    Velocity2D VelocityField::sample_at_horizontal_face(int i, int j) const {
        /// Get velocity vector at the horizontal face on the bottom of cell [i, j]
        
        /// interpolates either the 2 or 4 closes horizontal faces to obtain v
        double u = (read_u_or_outside(i, j) + read_u_or_outside(i + 1, j) +
            read_u_or_outside(i, j - 1) + read_u_or_outside(i + 1, j - 1)) / 4.0;
        double v = read_v_or_outside(i, j);

        return {u, v};
    }

    Velocity2D VelocityField::sample_at_corner(int i, int j) const {
        /// returns the velocity vector at the bottom left corenr of cell [i, j]

        Velocity2D up    = sample_at_vertical_face(i, j);
        Velocity2D right = sample_at_horizontal_face(i, j);
        Velocity2D down  = sample_at_vertical_face(i, j - 1);
        Velocity2D left  = sample_at_horizontal_face(i - 1, j);

        double u = (up.u + down.u) / 2.0;
        double v = (left.v + right.v) / 2.0;

        return {u, v};
    }

    Velocity2D VelocityField::sample_at_coordinates(double x, double y) const {
        /// uses interpolation of closest 4 walls 
        /// to return the velocity vector at given coordinates

        try {
            check_coordinates_bounds(x, y);
        }
        catch(const std::out_of_range& e) {
            return outside_velocity_;
        }

        int i = std::floor(x / resolution_);
        int j = std::floor(y / resolution_);
        
        double alpha_x = (x - static_cast<double>(i) * resolution_) / resolution_; 
        double alpha_y = (y - static_cast<double>(j) * resolution_) / resolution_; 

        Velocity2D ld = sample_at_corner(i, j);
        Velocity2D lu = sample_at_corner(i, j + 1);
        Velocity2D ru = sample_at_corner(i + 1, j + 1);
        Velocity2D rd = sample_at_corner(i + 1, j);

        double u = (1 - alpha_x) * (1 - alpha_y) * ld.u + 
                    alpha_x * (1 - alpha_y) * rd.u +
                    (1 - alpha_x) * alpha_y * lu.u +
                    alpha_x * alpha_y * ru.u;
        
        double v = (1 - alpha_x) * (1 - alpha_y) * ld.v + 
                    alpha_x * (1 - alpha_y) * rd.v +
                    (1 - alpha_x) * alpha_y * lu.v +
                    alpha_x * alpha_y * ru.v;

        return {u, v};
    }

}