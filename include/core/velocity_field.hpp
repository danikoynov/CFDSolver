#pragma once
#include <vector>

namespace cfd {

    constexpr double EPS = 1e-12;

    struct Velocity2D {
        double u;
        double v;
    };

    class VelocityField {
        /** 
         * Index pair (i, j) means the following things:
         *  - Cell at width i and height j
         *  - The left face of cell (i, j)
         *  - The bottom face of cell(i, j)
        */
        
        private:
            std::size_t width_, height_;
            double resolution_;
            Velocity2D outside_velocity_;
            std::vector<double> u_; // u_ is horizontal velocity field [width + 1, height]
            std::vector<double> v_; // v_ is vertical velocity field [width, height + 1]
        
            void check_u_bounds(int i, int j) const;
            void check_v_bounds(int i, int j) const;
            void check_coordinates_bounds(double x, double y) const;
        
            public:

            VelocityField(std::size_t width, std::size_t height, double resolution);
            
            double& get_u(int i, int j);
            const double& get_u(int i, int j) const;
            double& get_v(int i, int j);
            const double& get_v(int i, int j) const;

            double read_v_or_outside(int i, int j) const;
            double read_u_or_outside(int i, int j) const;

            Velocity2D sample_at_horizontal_face(int i, int j) const;
            Velocity2D sample_at_vertical_face(int i, int j) const;
            Velocity2D sample_at_coordinates(double x, double y) const;
            Velocity2D sample_at_corner(int i, int j) const;

    };
}