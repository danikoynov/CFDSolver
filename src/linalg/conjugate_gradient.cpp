#include "linalg/conjugate_gradient.hpp"
#include <cmath>

namespace cfd::linalg {

    Vector conjugate_gradient(const LinearOperator& A, const Vector& b) {

        Vector x(b.n());
        Vector r = b - A.apply(x);
        Vector p = r;

        const double TOL = 1e-9;
        
        auto within_tolerance = [TOL] (const Vector& v) {
            for (std::size_t i = 0; i < v.n(); i ++) {
                if (std::abs(v(i)) > TOL) {
                    return false;
                }
            }
            return true;
        };

        while(!within_tolerance(r)) {
            Vector Ap = A.apply(p);

            double rr = Vector::dot(r, r);
            double alpha = rr / Vector::dot(p, Ap);

            x = x + alpha * p;
            Vector r_new = r - alpha * Ap;

            double beta = Vector::dot(r_new, r_new) / rr;
            p = r_new + beta * p;
            r = r_new;
        }

        return x;
    }
    
}