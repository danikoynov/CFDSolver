#pragma once
#include <vector>

namespace cfd {

    class PressureField {

        private:
            std::size_t width_, height_;     
            std::vector<double> p_;
            void check_bounds(int i, int j) const;
        
            public:

            PressureField(std::size_t width, std::size_t height);

            double& get_p(int i, int j);
            const double& get_p(int i, int j) const;

    };
}