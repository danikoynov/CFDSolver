#include "linalg/gaussian_elimination.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using cfd::linalg::Matrix;
using cfd::linalg::Vector;
using cfd::linalg::gaussian_elimination;

namespace {

    constexpr double TOL = 1e-9;

    bool approx_equal(double a, double b, double tol = TOL) {
        return std::abs(a - b) < tol;
    }

    void require(bool condition, const std::string& message) {
        if (!condition) {
            throw std::runtime_error(message);
        }
    }

    void require_vector_close(
        const Vector& v,
        std::initializer_list<double> expected,
        const std::string& test_name
    ) {
        require(v.n() == expected.size(), test_name + ": wrong vector size");

        std::size_t i = 0;
        for (double x : expected) {
            if (!approx_equal(v(i), x)) {
                throw std::runtime_error(
                    test_name + ": wrong value at index " + std::to_string(i)
                );
            }
            ++i;
        }
    }

    void test_dimension_mismatch_throws() {
        Matrix A(2, 2);
        Vector b(3);

        A(0, 0) = 1.0; A(0, 1) = 2.0;
        A(1, 0) = 3.0; A(1, 1) = 4.0;

        b(0) = 5.0; b(1) = 6.0; b(2) = 7.0;

        bool threw = false;
        try {
            gaussian_elimination(A, b);
        } catch (const std::runtime_error&) {
            threw = true;
        }

        require(threw, "test_dimension_mismatch_throws failed");
    }

    void test_more_unknowns_than_equations_throws() {
        Matrix A(2, 3);
        Vector b(2);

        A(0, 0) = 1.0; A(0, 1) = 2.0; A(0, 2) = 3.0;
        A(1, 0) = 4.0; A(1, 1) = 5.0; A(1, 2) = 6.0;

        b(0) = 7.0; b(1) = 8.0;

        bool threw = false;
        try {
            gaussian_elimination(A, b);
        } catch (const std::runtime_error&) {
            threw = true;
        }

        require(threw, "test_more_unknowns_than_equations_throws failed");
    }

    void test_missing_pivot_throws() {
        // Dependent system:
        // x + y = 2
        // 2x + 2y = 4
        Matrix A(2, 2);
        Vector b(2);

        A(0, 0) = 1.0; A(0, 1) = 1.0;
        A(1, 0) = 2.0; A(1, 1) = 2.0;

        b(0) = 2.0;
        b(1) = 4.0;

        bool threw = false;
        try {
            gaussian_elimination(A, b);
        } catch (const std::runtime_error&) {
            threw = true;
        }

        require(threw, "test_missing_pivot_throws failed");
    }

    void test_no_solution_throws() {
        // Overdetermined inconsistent system:
        // x + y = 2
        // x - y = 0
        // 2x + 2y = 5   <- inconsistent with first equation
        Matrix A(3, 2);
        Vector b(3);

        A(0, 0) = 1.0; A(0, 1) = 1.0;
        A(1, 0) = 1.0; A(1, 1) = -1.0;
        A(2, 0) = 2.0; A(2, 1) = 2.0;

        b(0) = 2.0;
        b(1) = 0.0;
        b(2) = 5.0;

        bool threw = false;
        try {
            gaussian_elimination(A, b);
        } catch (const std::runtime_error&) {
            threw = true;
        }

        require(threw, "test_no_solution_throws failed");
    }

    void test_1x1_system() {
        // 4x = 20
        Matrix A(1, 1);
        Vector b(1);

        A(0, 0) = 4.0;
        b(0) = 20.0;

        Vector x = gaussian_elimination(A, b);
        require_vector_close(x, {5.0}, "test_1x1_system");
    }

    void test_2x2_system() {
        // x + y = 3
        // 2x + y = 4
        // solution: x=1, y=2
        Matrix A(2, 2);
        Vector b(2);

        A(0, 0) = 1.0; A(0, 1) = 1.0;
        A(1, 0) = 2.0; A(1, 1) = 1.0;

        b(0) = 3.0;
        b(1) = 4.0;

        Vector x = gaussian_elimination(A, b);
        require_vector_close(x, {1.0, 2.0}, "test_2x2_system");
    }

    void test_3x3_system() {
        // x + y + z = 6
        // 2x - y + z = 3
        // x + 2y - z = 3
        // solution: x=9/7, y=15/7, z=18/7
        Matrix A(3, 3);
        Vector b(3);

        A(0, 0) = 1.0; A(0, 1) = 1.0;  A(0, 2) = 1.0;
        A(1, 0) = 2.0; A(1, 1) = -1.0; A(1, 2) = 1.0;
        A(2, 0) = 1.0; A(2, 1) = 2.0;  A(2, 2) = -1.0;

        b(0) = 6.0;
        b(1) = 3.0;
        b(2) = 3.0;

        Vector x = gaussian_elimination(A, b);
        require_vector_close(
            x,
            {9.0 / 7.0, 15.0 / 7.0, 18.0 / 7.0},
            "test_3x3_system"
        );
    }

    void test_row_swap_needed() {
        // 0x + y = 2
        // x + y = 3
        // solution: x=1, y=2
        Matrix A(2, 2);
        Vector b(2);

        A(0, 0) = 0.0; A(0, 1) = 1.0;
        A(1, 0) = 1.0; A(1, 1) = 1.0;

        b(0) = 2.0;
        b(1) = 3.0;

        Vector x = gaussian_elimination(A, b);
        require_vector_close(x, {1.0, 2.0}, "test_row_swap_needed");
    }

    void test_overdetermined_consistent_system() {
        // x + y = 2
        // x - y = 0
        // 2x = 2
        // solution: x=1, y=1
        Matrix A(3, 2);
        Vector b(3);

        A(0, 0) = 1.0; A(0, 1) = 1.0;
        A(1, 0) = 1.0; A(1, 1) = -1.0;
        A(2, 0) = 2.0; A(2, 1) = 0.0;

        b(0) = 2.0;
        b(1) = 0.0;
        b(2) = 2.0;

        Vector x = gaussian_elimination(A, b);
        require_vector_close(x, {1.0, 1.0}, "test_overdetermined_consistent_system");
    }

    void test_decimal_coefficients() {
        // 0.5x + 1.5y = 5
        // 1.0x - 1.0y = -1
        // solution: x=1, y=2
        Matrix A(2, 2);
        Vector b(2);

        A(0, 0) = 0.5; A(0, 1) = 1.5;
        A(1, 0) = 1.0; A(1, 1) = -1.0;

        b(0) = 3.5;
        b(1) = -1.0;

        Vector x = gaussian_elimination(A, b);
        require_vector_close(x, {1.0, 2.0}, "test_decimal_coefficients");
    }

} // namespace

int main() {
    int passed = 0;
    int failed = 0;

    auto run = [&](void (*test)(), const char* name) {
        try {
            test();
            std::cout << "[PASS] " << name << '\n';
            ++passed;
        } catch (const std::exception& e) {
            std::cout << "[FAIL] " << name << " -> " << e.what() << '\n';
            ++failed;
        }
    };

    run(test_dimension_mismatch_throws, "dimension mismatch throws");
    run(test_more_unknowns_than_equations_throws, "more unknowns than equations throws");
    run(test_missing_pivot_throws, "missing pivot throws");
    run(test_no_solution_throws, "no solution throws");

    run(test_1x1_system, "1x1 system");
    run(test_2x2_system, "2x2 system");
    run(test_3x3_system, "3x3 system");
    run(test_row_swap_needed, "row swap needed");
    run(test_overdetermined_consistent_system, "overdetermined consistent system");
    run(test_decimal_coefficients, "decimal coefficients");

    std::cout << "\nPassed: " << passed << '\n';
    std::cout << "Failed: " << failed << '\n';

    return failed == 0 ? 0 : 1;
}