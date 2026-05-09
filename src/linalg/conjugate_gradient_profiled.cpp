#include "linalg/conjugate_gradient.hpp"
#include <cmath>
#include <chrono>
#include <iostream>
#include <iomanip>


namespace cfd::linalg {
 
    Vector conjugate_gradient_profiled(const LinearOperator& A, const Vector& b) {
        using clock = std::chrono::steady_clock;

        auto elapsed_ms = [](const clock::time_point& a, const clock::time_point& b) {
            return std::chrono::duration<double, std::milli>(b - a).count();
        };

        auto measure = [&](double& total_ms, int& calls, auto&& fn) {
            auto t0 = clock::now();
            fn();
            auto t1 = clock::now();
            total_ms += elapsed_ms(t0, t1);
            ++calls;
        };

        double total_ms = 0.0;
        double init_ms = 0.0;
        double tol_check_ms = 0.0;
        double apply_ms = 0.0;
        double dot_ms = 0.0;
        double x_update_ms = 0.0;
        double r_update_ms = 0.0;
        double p_update_ms = 0.0;

        int init_calls = 0;
        int tol_check_calls = 0;
        int apply_calls = 0;
        int dot_calls = 0;
        int x_update_calls = 0;
        int r_update_calls = 0;
        int p_update_calls = 0;

        const double TOL = 1e-3;

        auto within_tolerance = [TOL](const Vector& v) {
            for (std::size_t i = 0; i < v.n(); i++) {
                if (std::abs(v(i)) > TOL) {
                    return false;
                }
            }
            return true;
        };

        auto total_start = clock::now();

        Vector x(b.n());
        Vector r(b.n());
        Vector p(b.n());

        measure(init_ms, init_calls, [&] {
            Vector Ax0(b.n());

            measure(apply_ms, apply_calls, [&] {
                Ax0 = A.apply(x);
            });

            r = b - Ax0;
            p = r;
        });

        int iterations = 0;

        while (true) {
            bool done = false;
            measure(tol_check_ms, tol_check_calls, [&] {
                done = within_tolerance(r);
            });
            if (done) {
                break;
            }

            ++iterations;

            Vector Ap(b.n());
            measure(apply_ms, apply_calls, [&] {
                Ap = A.apply(p);
            });

            double rr = 0.0;
            double pAp = 0.0;

            measure(dot_ms, dot_calls, [&] {
                rr = Vector::dot(r, r);
            });

            measure(dot_ms, dot_calls, [&] {
                pAp = Vector::dot(p, Ap);
            });

            double alpha = rr / pAp;

            measure(x_update_ms, x_update_calls, [&] {
                x += alpha * p;
            });

            Vector r_new(b.n());
            measure(r_update_ms, r_update_calls, [&] {
                r_new = r - alpha * Ap;
            });

            double rnew_rnew = 0.0;
            measure(dot_ms, dot_calls, [&] {
                rnew_rnew = Vector::dot(r_new, r_new);
            });

            double beta = rnew_rnew / rr;

            measure(p_update_ms, p_update_calls, [&] {
                p = r_new + beta * p;
            });

            r = std::move(r_new);
        }

        auto total_end = clock::now();
        total_ms = elapsed_ms(total_start, total_end);

        auto print_stat = [&](const char* name, double ms, int calls) {
            double pct = (total_ms > 0.0) ? (100.0 * ms / total_ms) : 0.0;

            std::cerr
                << std::left << std::setw(22) << name
                << " | total: " << std::setw(10) << ms << " ms"
                << " | pct: " << std::setw(8) << pct << " %"
                << " | calls: " << std::setw(6) << calls;

            if (calls > 0) {
                std::cerr << " | avg/call: " << (ms / calls) << " ms";
            }

            std::cerr << '\n';
        };

        std::cerr << std::fixed << std::setprecision(3);
        std::cerr << "\n=== CG Profile ===\n";
        std::cerr << "Iterations: " << iterations << '\n';
        std::cerr << "Total CG time: " << total_ms << " ms\n";
        print_stat("initialization", init_ms, init_calls);
        print_stat("tolerance check", tol_check_ms, tol_check_calls);
        print_stat("A.apply", apply_ms, apply_calls);
        print_stat("dot products", dot_ms, dot_calls);
        print_stat("x update", x_update_ms, x_update_calls);
        print_stat("r update", r_update_ms, r_update_calls);
        print_stat("p update", p_update_ms, p_update_calls);
        std::cerr << "==================\n";


        return x;
    }
        
}