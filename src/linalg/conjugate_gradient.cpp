#include "linalg/conjugate_gradient.hpp"

#include <cmath>
#include <stdexcept>

namespace cfd::linalg {
 
    Vector conjugate_gradient(const PoissonOperator& A, const Vector& b) {
        const double TOL = 1e-3;
        const int MAX_ITERATIONS = 10000;

        auto within_tolerance = [TOL](const Vector& v) {
            for (std::size_t i = 0; i < v.n(); i++) {
                if (std::abs(v(i)) > TOL) {
                    return false;
                }
            }
            return true;
        };

        Vector x(b.n());

        Vector Ax = A.apply(x);
        Vector r = b - Ax;

        Vector z = A.apply_preconditioner(r);
        Vector p = z;

        double rz = Vector::dot(r, z);

        int iterations = 0;

        while (!within_tolerance(r)) {
            if (iterations >= MAX_ITERATIONS) {
                throw std::runtime_error("Conjugate gradient failed to converge.");
            }

            iterations++;

            Vector Ap = A.apply(p);

            const double pAp = Vector::dot(p, Ap);

            if (std::abs(pAp) < 1e-20) {
                throw std::runtime_error("Conjugate gradient breakdown: pAp is too small.");
            }

            const double alpha = rz / pAp;

            for (std::size_t i = 0; i < b.n(); i++) {
                x(i) += alpha * p(i);
                r(i) -= alpha * Ap(i);
            }

            Vector z_new = A.apply_preconditioner(r);

            const double rz_new = Vector::dot(r, z_new);

            if (std::abs(rz) < 1e-20) {
                throw std::runtime_error("Conjugate gradient breakdown: rz is too small.");
            }

            const double beta = rz_new / rz;

            for (std::size_t i = 0; i < b.n(); i++) {
                p(i) = z_new(i) + beta * p(i);
            }

            z = std::move(z_new);
            rz = rz_new;
        }

        return x;
    }
        
}