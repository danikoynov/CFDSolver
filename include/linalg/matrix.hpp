#pragma once
#include<vector>

namespace cfd::linalg {

    class Matrix{
        
        private:
            std::size_t rows_, cols_;
            std::vector<double> data_; /// Matrix stored in one vector
                                       /// for less allocation overhead
        
        public:
            Matrix(std::size_t rows, std::size_t cols);
            double& operator()(std::size_t i, std::size_t j);
            const double& operator()(std::size_t i, std::size_t j) const;
            void swap_rows(std::size_t i, std::size_t j);
            std::size_t rows() const;
            std::size_t cols() const;
    };   
}