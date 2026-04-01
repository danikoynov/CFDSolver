#include<vector>

namespace cfd::linalg {

    class Matrix{
        
        private:
            std::size_t rows_, cols_;
            std::vector<double> data_; /// Matrix stored in one vector
                                       /// for less allocation overhead
        
        public:
            Matrix(std::size_t rows, std::size_t cols);
            double& operator()(size_t i, size_t j);

    }   
}