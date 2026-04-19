#include "linalg/conjugate_gradient.hpp"
#include "linalg/linear_operator.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using cfd::linalg::LinearOperator;
using cfd::linalg::Matrix;
using cfd::linalg::Vector;
using cfd::linalg::conjugate_gradient;

namespace {

    constexpr double TOL = 1e-6;

    bool approx_equal(double a, double b, double tol = TOL) {
        return std::abs(a - b) < tol;
    }

    void require(bool condition, const std::string& message) {
        if (!condition) {
            throw std::runtime_error(message);
        }
    }

    LinearOperator make_operator(const Matrix& A) {
        LinearOperator op(A.cols(), A.rows());

        for (std::size_t i = 0; i < A.rows(); i++) {
            for (std::size_t j = 0; j < A.cols(); j++) {
                if (std::abs(A(i, j)) > 0.0) {
                    op(static_cast<int>(i), static_cast<int>(j)) = A(i, j);
                }
            }
        }

        op.init_csr();
        return op;
    }

    Vector matvec(const Matrix& A, const Vector& x) {
        require(A.cols() == x.n(), "matvec: dimension mismatch");

        Vector res(A.rows());

        for (std::size_t i = 0; i < A.rows(); i++) {
            for (std::size_t j = 0; j < A.cols(); j++) {
                res(i) += A(i, j) * x(j);
            }
        }

        return res;
    }

    void require_vector_close(
        const Vector& v,
        std::initializer_list<double> expected,
        const std::string& test_name,
        double tol = TOL
    ) {
        require(v.n() == expected.size(), test_name + ": wrong vector size");

        std::size_t i = 0;
        for (double x : expected) {
            if (!approx_equal(v(i), x, tol)) {
                throw std::runtime_error(
                    test_name + ": wrong value at index " + std::to_string(i)
                );
            }
            ++i;
        }
    }

    void require_residual_small(
        const Matrix& A,
        const Vector& x,
        const Vector& b,
        const std::string& test_name,
        double tol = TOL
    ) {
        Vector Ax = matvec(A, x);

        require(Ax.n() == b.n(), test_name + ": residual size mismatch");

        for (std::size_t i = 0; i < b.n(); i++) {
            if (std::abs(Ax(i) - b(i)) > tol) {
                throw std::runtime_error(
                    test_name + ": residual too large at index " + std::to_string(i)
                );
            }
        }
    }

    void test_dimension_mismatch_throws() {
        Matrix A(2, 2);
        Vector b(3);

        A(0, 0) = 4.0; A(0, 1) = 1.0;
        A(1, 0) = 1.0; A(1, 1) = 3.0;

        b(0) = 1.0; b(1) = 2.0; b(2) = 3.0;

        LinearOperator op = make_operator(A);

        bool threw = false;
        try {
            conjugate_gradient(op, b);
        } catch (const std::runtime_error&) {
            threw = true;
        }

        require(threw, "test_dimension_mismatch_throws failed");
    }

    void test_non_square_matrix_throws() {
        Matrix A(2, 3);
        Vector b(2);

        A(0, 0) = 1.0; A(0, 1) = 0.0; A(0, 2) = 0.0;
        A(1, 0) = 0.0; A(1, 1) = 1.0; A(1, 2) = 0.0;

        b(0) = 1.0; b(1) = 2.0;

        LinearOperator op = make_operator(A);

        bool threw = false;
        try {
            conjugate_gradient(op, b);
        } catch (const std::runtime_error&) {
            threw = true;
        }

        require(threw, "test_non_square_matrix_throws failed");
    }

    void test_zero_rhs_returns_zero() {
        Matrix A(3, 3);
        Vector b(3);

        A(0, 0) = 4.0; A(0, 1) = 1.0; A(0, 2) = 0.0;
        A(1, 0) = 1.0; A(1, 1) = 3.0; A(1, 2) = 1.0;
        A(2, 0) = 0.0; A(2, 1) = 1.0; A(2, 2) = 2.0;

        b(0) = 0.0; b(1) = 0.0; b(2) = 0.0;

        LinearOperator op = make_operator(A);
        Vector x = conjugate_gradient(op, b);

        require_vector_close(x, {0.0, 0.0, 0.0}, "test_zero_rhs_returns_zero");
    }

    void test_1x1_system() {
        Matrix A(1, 1);
        Vector b(1);

        A(0, 0) = 4.0;
        b(0) = 20.0;

        LinearOperator op = make_operator(A);
        Vector x = conjugate_gradient(op, b);

        require_vector_close(x, {5.0}, "test_1x1_system");
    }

    void test_2x2_spd_system() {
        Matrix A(2, 2);
        Vector b(2);

        A(0, 0) = 4.0; A(0, 1) = 1.0;
        A(1, 0) = 1.0; A(1, 1) = 3.0;

        b(0) = 6.0;
        b(1) = 7.0;

        LinearOperator op = make_operator(A);
        Vector x = conjugate_gradient(op, b);

        require_vector_close(x, {1.0, 2.0}, "test_2x2_spd_system");
        require_residual_small(A, x, b, "test_2x2_spd_system");
    }

    void test_3x3_spd_system() {
        Matrix A(3, 3);
        Vector b(3);

        A(0, 0) = 4.0; A(0, 1) = 1.0; A(0, 2) = 0.0;
        A(1, 0) = 1.0; A(1, 1) = 3.0; A(1, 2) = 1.0;
        A(2, 0) = 0.0; A(2, 1) = 1.0; A(2, 2) = 2.0;

        b(0) = 6.0;
        b(1) = 10.0;
        b(2) = 8.0;

        LinearOperator op = make_operator(A);
        Vector x = conjugate_gradient(op, b);

        require_vector_close(x, {1.0, 2.0, 3.0}, "test_3x3_spd_system");
        require_residual_small(A, x, b, "test_3x3_spd_system");
    }

    void test_diagonal_spd_system() {
        Matrix A(4, 4);
        Vector b(4);

        A(0, 0) = 2.0; A(0, 1) = 0.0; A(0, 2) = 0.0; A(0, 3) = 0.0;
        A(1, 0) = 0.0; A(1, 1) = 3.0; A(1, 2) = 0.0; A(1, 3) = 0.0;
        A(2, 0) = 0.0; A(2, 1) = 0.0; A(2, 2) = 4.0; A(2, 3) = 0.0;
        A(3, 0) = 0.0; A(3, 1) = 0.0; A(3, 2) = 0.0; A(3, 3) = 5.0;

        b(0) = 2.0;
        b(1) = 6.0;
        b(2) = 12.0;
        b(3) = 20.0;

        LinearOperator op = make_operator(A);
        Vector x = conjugate_gradient(op, b);

        require_vector_close(x, {1.0, 2.0, 3.0, 4.0}, "test_diagonal_spd_system");
        require_residual_small(A, x, b, "test_diagonal_spd_system");
    }

    void test_poisson_like_tridiagonal_system() {
        Matrix A(4, 4);
        Vector b(4);

        A(0, 0) = 2.0;  A(0, 1) = -1.0; A(0, 2) = 0.0;  A(0, 3) = 0.0;
        A(1, 0) = -1.0; A(1, 1) = 2.0;  A(1, 2) = -1.0; A(1, 3) = 0.0;
        A(2, 0) = 0.0;  A(2, 1) = -1.0; A(2, 2) = 2.0;  A(2, 3) = -1.0;
        A(3, 0) = 0.0;  A(3, 1) = 0.0;  A(3, 2) = -1.0; A(3, 3) = 2.0;

        b(0) = 0.0;
        b(1) = 0.0;
        b(2) = 0.0;
        b(3) = 5.0;

        LinearOperator op = make_operator(A);
        Vector x = conjugate_gradient(op, b);

        require_vector_close(
            x,
            {1.0, 2.0, 3.0, 4.0},
            "test_poisson_like_tridiagonal_system"
        );
        require_residual_small(A, x, b, "test_poisson_like_tridiagonal_system");
    }

    void test_decimal_spd_system() {
        Matrix A(2, 2);
        Vector b(2);

        A(0, 0) = 1.5; A(0, 1) = 0.5;
        A(1, 0) = 0.5; A(1, 1) = 2.0;

        b(0) = 2.5;
        b(1) = 4.5;

        LinearOperator op = make_operator(A);
        Vector x = conjugate_gradient(op, b);

        require_vector_close(x, {1.0, 2.0}, "test_decimal_spd_system");
        require_residual_small(A, x, b, "test_decimal_spd_system");
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
    run(test_non_square_matrix_throws, "non-square matrix throws");

    run(test_zero_rhs_returns_zero, "zero rhs returns zero");
    run(test_1x1_system, "1x1 system");
    run(test_2x2_spd_system, "2x2 spd system");
    run(test_3x3_spd_system, "3x3 spd system");
    run(test_diagonal_spd_system, "diagonal spd system");
    run(test_poisson_like_tridiagonal_system, "poisson-like tridiagonal system");
    run(test_decimal_spd_system, "decimal spd system");

    std::cout << "\nPassed: " << passed << '\n';
    std::cout << "Failed: " << failed << '\n';

    return failed == 0 ? 0 : 1;
}