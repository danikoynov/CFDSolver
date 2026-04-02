#pragma once
#include<vector>

namespace cfd::linalg {

    class Vector{
        
        private:
            std::size_t n_;
            std::vector<double> data_;
        
        public:
            Vector(std::size_t n);
            double& operator()(std::size_t i);
            const double& operator()(std::size_t i) const;
            std::size_t n() const;

    };   
}