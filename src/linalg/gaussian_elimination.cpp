#include "linalg/gaussian_elimination.hpp"
#include <stdexcept>
#include <utility>
#include <cmath>

namespace cfd::linalg {

    constexpr double EPS = 1e-12;

    Vector gaussian_elimination(Matrix A, Vector b) {
        /// Solves a system of linaer equations in the form Ax=b

        if (A.rows() != b.n()) {
            throw std::runtime_error(
                "Error during Gaussian elimination:\n"
                "number of equations does not match "
                "number of reults of equations."
            );
        }

        if (A.cols() > A.rows()) {
            throw std::runtime_error(
                "Error during Gaussian elimination:\n"
                "number of unknowns is more than number of equations."
            );
        }

        /// Perform row reduction
        for (std::size_t i = 0; i < A.cols(); i ++) {
            
            std::size_t pivot = i + 1;

            while((std::abs(A(i, i)) < EPS) && pivot < A.rows()) {
                if (!(std::abs(A(pivot, i)) < EPS)) {
                    A.swap_rows(i, pivot);
                    std::swap(b(i), b(pivot));
                    break;
                }
                pivot ++;
            }
            
            if (std::abs(A(i, i)) < EPS) {
                throw std::runtime_error(
                    "Error during Gaussian elimination: \n"
                    "system of equations has a missing pivot."
                );
            }

            b(i) /= A(i, i);
            for (std::size_t j = A.cols(); j-- > i;) { /// Normalize row
                A(i, j) /= A(i, i);
            }
            
            
            for (std::size_t k = i + 1; k < A.rows(); k ++) {
                b(k) -= A(k, i) * b(i);
                for (std::size_t j = A.cols(); j-- > i;) 
                    A(k, j) -= A(k, i) * A(i, j);
            }
                    
        }

        for (std::size_t i = A.cols(); i < A.rows(); i ++)
            if (!(std::abs(b(i)) < EPS))
                throw std::runtime_error(
                    "Error during Gaussian elimination: \n"
                    "system of equations has no solution."
                );


        Vector sol(A.cols());

        for (std::size_t i = A.cols() ; i-- > 0;) {
            double value = b(i);
            
            for (std::size_t j = A.cols() - 1; j > i; j --) {
                value -= sol(j) * A(i, j); 
            }

            sol(i) = value;
        }

        return sol;
    }
}