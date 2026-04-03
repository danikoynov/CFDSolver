#pragma once
#include <vector>

namespace cfd {

    class PressureField {

        private:
            std::size_t width_, height_;
            const double outside_pressure_;     
            std::vector<double> p_;
            void check_bounds(int i, int j) const;
        
            public:

            PressureField(std::size_t width, std::size_t height, 
                double outside_pressure);

            double& get_p(int i, int j);
            const double& get_p(int i, int j) const;

            double read_p_or_outside(int i, int j) const;

    };
}