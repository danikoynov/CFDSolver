#pragma once
#include "linalg/matrix.hpp"
#include "linalg/vector.hpp"
#include <vector>

namespace cfd::linalg {
    
    inline constexpr double EPS = 1e-9;

    class PoissonOperator {

        public:
            PoissonOperator(
                std::size_t width,
                std::size_t height);
                
            Vector apply(const Vector& v) const;
            // Apply the operator on vector v

            Vector apply_preconditioner(const Vector& r) const;
            // Apply the IC preconditioner on vector r

            double& operator()(int x, int y);
            // Provide mutable access to entry (x, y)

            const double& operator()(int x, int y) const;
            // Provide immutable access to entry (x, y)

            double get(int i, int j) const; 
            // Return the value of entry (i, j) or 0.0 if not present

            double get_diag(int k) const;   // Return the value of entry (k, k)
            double get_bottom(int k) const; // Return the value of entry (k, k - width)
            double get_left(int k) const;   // Return the value of entry (k, k - 1)
            double get_top(int k) const;    // Return the value of entry (k, k + width)
            double get_right(int k) const;  // Return the value of entry (k, k + 1)
            
            void finalize(); // Run the sterilizer and initializers

        private:

            void init_preconditioner(); // Initialize Incomplete Cholesky preconditioner
            void init_csr();            // Initialize compressed sparse row representation
            void sterilize();           // Remove zero entries in the operator

            std::size_t width_;  // Width of grid
            std::size_t height_; // Height of grid
            std::size_t size_;   // Size of the Poisson operator (square matrix)

            std::vector<std::vector<std::pair<std::size_t, double>>> assembly_rows_; 
            // Assembly storage: each row contains (col, value) pairs

            std::vector<std::pair<std::size_t, double>> csr_entries_; 
            // CSR storage: flattened (column, value) entries used during apply()

            std::vector<std::size_t> row_ptr_; 
            // CSR row pointers: row i occupies entries inside [row_ptr_[i], row_ptr_[i + 1])
            
            std::vector<double> ic_inv_diag_; // Incomplete Cholesky inverse diagonal:  1 / L[k, k]
            std::vector<double> ic_bottom_;   // Incomplete Cholesky lower bottom entry: L[k, k - width]
            std::vector<double> ic_left_;     // Incomplete Cholesky lower left entry:  L[k, k - 1]
    };

}