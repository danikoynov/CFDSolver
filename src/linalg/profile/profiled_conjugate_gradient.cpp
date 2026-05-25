#include "linalg/conjugate_gradient.hpp"

#include <cmath>
#include <stdexcept>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <filesystem>
#include <optional>

namespace cfd::linalg {

    Vector profiled_conjugate_gradient(
        const PoissonOperator& A,
        const Vector& b,
        bool save_data,
        std::optional<std::filesystem::path> output_dir
    ) {
        if (save_data && !output_dir.has_value()) {
            throw std::invalid_argument(
                "Cannot save CG profiling data without providing an output directory."
            );
        }

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

        auto avg_ms = [](double total_ms, int calls) {
            if (calls == 0) {
                return 0.0;
            }

            return total_ms / static_cast<double>(calls);
        };

        double total_ms = 0.0;
        double init_ms = 0.0;
        double apply_precond_ms = 0.0;
        double dot_ms = 0.0;

        int apply_precond_calls = 0;
        int dot_calls = 0;

        auto start = clock::now();

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

        auto start_init = clock::now();

        Vector x(b.n());

        Vector Ax = A.apply(x);
        Vector r = b - Ax;

        Vector z = A.apply_preconditioner(r);
        Vector p = z;

        double rz = Vector::dot(r, z);

        int iterations = 0;

        auto end_init = clock::now();
        init_ms = elapsed_ms(start_init, end_init);

        while (!within_tolerance(r)) {
            if (iterations >= MAX_ITERATIONS) {
                throw std::runtime_error("Conjugate gradient failed to converge.");
            }

            iterations++;

            Vector Ap = A.apply(p);

            double pAp = 0.0;

            measure(dot_ms, dot_calls, [&]() {
                pAp = Vector::dot(p, Ap);
            });

            if (std::abs(pAp) < 1e-20) {
                throw std::runtime_error("Conjugate gradient breakdown: pAp is too small.");
            }

            const double alpha = rz / pAp;

            for (std::size_t i = 0; i < b.n(); i++) {
                x(i) += alpha * p(i);
                r(i) -= alpha * Ap(i);
            }

            Vector z_new(b.n());

            measure(apply_precond_ms, apply_precond_calls, [&]() {
                z_new = A.apply_preconditioner(r);
            });

            double rz_new = 0.0;

            measure(dot_ms, dot_calls, [&]() {
                rz_new = Vector::dot(r, z_new);
            });

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

        auto end = clock::now();
        total_ms = elapsed_ms(start, end);

        auto print_measurement = [&](
            const std::string& name,
            double total_ms,
            int calls
        ) {
            std::cout << name << " | "
                      << "total: " << std::setw(7) << total_ms << " ms | "
                      << "calls: " << std::setw(3) << calls << " | "
                      << "avg/call " << avg_ms(total_ms, calls) << " ms"
                      << std::endl;
        };

        std::cout << "=== Conjugate Gradient Profile ===" << std::endl;
        std::cout << std::fixed << std::setprecision(3);
        std::cout << "Total Time: " << std::setw(8) << total_ms << " ms" << std::endl;
        std::cout << "Iterations: " << std::setw(8) << iterations << std::endl;
        std::cout << "Init Time : " << std::setw(8) << init_ms << " ms" << std::endl;

        print_measurement("Precond", apply_precond_ms, apply_precond_calls);
        print_measurement("Dot    ", dot_ms, dot_calls);

        std::cout << std::endl;

        if (save_data) {
            namespace fs = std::filesystem;

            fs::create_directories(*output_dir);

            const fs::path csv_path = *output_dir / "dummy_cg_profile.csv";

            const bool write_header =
                !fs::exists(csv_path) || fs::file_size(csv_path) == 0;

            std::ofstream file(csv_path, std::ios::app);

            if (!file) {
                throw std::runtime_error(
                    "Could not open CG profiling CSV file: " + csv_path.string()
                );
            }

            if (write_header) {
                file << "iterations,"
                     << "total_ms,"
                     << "init_ms,"
                     << "preconditioner_total_ms,"
                     << "preconditioner_calls,"
                     << "preconditioner_avg_ms,"
                     << "dot_total_ms,"
                     << "dot_calls,"
                     << "dot_avg_ms\n";
            }

            file << iterations << ","
                 << total_ms << ","
                 << init_ms << ","
                 << apply_precond_ms << ","
                 << apply_precond_calls << ","
                 << avg_ms(apply_precond_ms, apply_precond_calls) << ","
                 << dot_ms << ","
                 << dot_calls << ","
                 << avg_ms(dot_ms, dot_calls)
                 << "\n";
        }

        return x;
    }

}