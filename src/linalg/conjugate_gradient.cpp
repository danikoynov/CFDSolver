#include "linalg/conjugate_gradient.hpp"
#include <cmath>

namespace cfd::linalg {

    Vector conjugate_gradient(const LinearOperator& A, const Vector& b) {
        const double TOL = 1e-6;
        const std::size_t MAX_ITER = 10 * b.n();

        Vector x(b.n());
        Vector r = b - A.apply(x);
        Vector p = r;
        double last_rr = Vector::dot(r, r);

        auto within_tolerance = [TOL](const Vector& v) {
            for (std::size_t i = 0; i < v.n(); i++) {
                if (std::abs(v(i)) > TOL) return false;
            }
            return true;
        };

        std::size_t iters = 0;
        while (iters++ < MAX_ITER) {
            if (iters % 10 == 0 && within_tolerance(r)) break;

            Vector Ap    = A.apply(p);
            double alpha = last_rr / Vector::dot(p, Ap);

            x.axpy( alpha, p);
            r.axpy(-alpha, Ap);

            double new_rr = Vector::dot(r, r);
            double beta   = new_rr / last_rr;
            last_rr       = new_rr;

            p *= beta;
            p.axpy(1.0, r);
        }

        return x;
    }

}