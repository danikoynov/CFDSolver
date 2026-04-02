#include "core/pressure_field.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

using cfd::PressureField;

namespace {

    void require(bool condition, const std::string& message) {
        if (!condition) {
            throw std::runtime_error(message);
        }
    }

    void test_constructor_initializes_to_zero() {
        PressureField field(3, 2, 101325.0);

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 2; ++j) {
                require(
                    field.get_p(i, j) == 0.0,
                    "test_constructor_initializes_to_zero failed"
                );
            }
        }
    }

    void test_get_p_write_and_read() {
        PressureField field(4, 3, 101325.0);

        field.get_p(0, 0) = 1.5;
        field.get_p(2, 1) = -3.25;
        field.get_p(3, 2) = 8.0;

        require(field.get_p(0, 0) == 1.5, "test_get_p_write_and_read failed at (0, 0)");
        require(field.get_p(2, 1) == -3.25, "test_get_p_write_and_read failed at (2, 1)");
        require(field.get_p(3, 2) == 8.0, "test_get_p_write_and_read failed at (3, 2)");
    }

    void test_const_get_p_reads_correct_value() {
        PressureField field(2, 2, 50.0);
        field.get_p(1, 1) = 7.75;

        const PressureField& const_field = field;

        require(
            const_field.get_p(1, 1) == 7.75,
            "test_const_get_p_reads_correct_value failed"
        );
    }

    void test_check_bounds_accepts_valid_indices() {
        PressureField field(3, 4, 0.0);

        field.check_bounds(0, 0);
        field.check_bounds(2, 3);
        field.check_bounds(1, 2);
    }

    void test_check_bounds_throws_on_negative_i() {
        PressureField field(3, 4, 0.0);

        bool threw = false;
        try {
            field.check_bounds(-1, 0);
        } catch (const std::out_of_range&) {
            threw = true;
        }

        require(threw, "test_check_bounds_throws_on_negative_i failed");
    }

    void test_check_bounds_throws_on_negative_j() {
        PressureField field(3, 4, 0.0);

        bool threw = false;
        try {
            field.check_bounds(0, -1);
        } catch (const std::out_of_range&) {
            threw = true;
        }

        require(threw, "test_check_bounds_throws_on_negative_j failed");
    }

    void test_check_bounds_throws_on_i_equal_width() {
        PressureField field(3, 4, 0.0);

        bool threw = false;
        try {
            field.check_bounds(3, 0);
        } catch (const std::out_of_range&) {
            threw = true;
        }

        require(threw, "test_check_bounds_throws_on_i_equal_width failed");
    }

    void test_check_bounds_throws_on_j_equal_height() {
        PressureField field(3, 4, 0.0);

        bool threw = false;
        try {
            field.check_bounds(0, 4);
        } catch (const std::out_of_range&) {
            threw = true;
        }

        require(threw, "test_check_bounds_throws_on_j_equal_height failed");
    }

    void test_get_p_throws_on_out_of_bounds_access() {
        PressureField field(2, 2, 0.0);

        bool threw = false;
        try {
            field.get_p(2, 1);
        } catch (const std::out_of_range&) {
            threw = true;
        }

        require(threw, "test_get_p_throws_on_out_of_bounds_access failed");
    }

    void test_read_p_or_outside_returns_inside_value() {
        PressureField field(3, 3, 123.0);
        field.get_p(1, 2) = 42.5;

        require(
            field.read_p_or_outside(1, 2) == 42.5,
            "test_read_p_or_outside_returns_inside_value failed"
        );
    }

    void test_read_p_or_outside_returns_outside_pressure_left() {
        PressureField field(3, 3, 123.0);

        require(
            field.read_p_or_outside(-1, 1) == 123.0,
            "test_read_p_or_outside_returns_outside_pressure_left failed"
        );
    }

    void test_read_p_or_outside_returns_outside_pressure_right() {
        PressureField field(3, 3, 123.0);

        require(
            field.read_p_or_outside(3, 1) == 123.0,
            "test_read_p_or_outside_returns_outside_pressure_right failed"
        );
    }

    void test_read_p_or_outside_returns_outside_pressure_bottom() {
        PressureField field(3, 3, 123.0);

        require(
            field.read_p_or_outside(1, -1) == 123.0,
            "test_read_p_or_outside_returns_outside_pressure_bottom failed"
        );
    }

    void test_read_p_or_outside_returns_outside_pressure_top() {
        PressureField field(3, 3, 123.0);

        require(
            field.read_p_or_outside(1, 3) == 123.0,
            "test_read_p_or_outside_returns_outside_pressure_top failed"
        );
    }

    void test_values_are_stored_independently() {
        PressureField field(3, 2, 0.0);

        field.get_p(0, 0) = 1.0;
        field.get_p(0, 1) = 2.0;
        field.get_p(1, 0) = 3.0;
        field.get_p(1, 1) = 4.0;
        field.get_p(2, 0) = 5.0;
        field.get_p(2, 1) = 6.0;

        require(field.get_p(0, 0) == 1.0, "test_values_are_stored_independently failed at (0, 0)");
        require(field.get_p(0, 1) == 2.0, "test_values_are_stored_independently failed at (0, 1)");
        require(field.get_p(1, 0) == 3.0, "test_values_are_stored_independently failed at (1, 0)");
        require(field.get_p(1, 1) == 4.0, "test_values_are_stored_independently failed at (1, 1)");
        require(field.get_p(2, 0) == 5.0, "test_values_are_stored_independently failed at (2, 0)");
        require(field.get_p(2, 1) == 6.0, "test_values_are_stored_independently failed at (2, 1)");
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

    run(test_constructor_initializes_to_zero, "constructor initializes to zero");
    run(test_get_p_write_and_read, "get_p write and read");
    run(test_const_get_p_reads_correct_value, "const get_p reads correct value");

    run(test_check_bounds_accepts_valid_indices, "check_bounds accepts valid indices");
    run(test_check_bounds_throws_on_negative_i, "check_bounds throws on negative i");
    run(test_check_bounds_throws_on_negative_j, "check_bounds throws on negative j");
    run(test_check_bounds_throws_on_i_equal_width, "check_bounds throws on i equal width");
    run(test_check_bounds_throws_on_j_equal_height, "check_bounds throws on j equal height");
    run(test_get_p_throws_on_out_of_bounds_access, "get_p throws on out of bounds access");

    run(test_read_p_or_outside_returns_inside_value, "read_p_or_outside returns inside value");
    run(test_read_p_or_outside_returns_outside_pressure_left, "read_p_or_outside returns outside pressure left");
    run(test_read_p_or_outside_returns_outside_pressure_right, "read_p_or_outside returns outside pressure right");
    run(test_read_p_or_outside_returns_outside_pressure_bottom, "read_p_or_outside returns outside pressure bottom");
    run(test_read_p_or_outside_returns_outside_pressure_top, "read_p_or_outside returns outside pressure top");

    run(test_values_are_stored_independently, "values are stored independently");

    std::cout << "\nPassed: " << passed << '\n';
    std::cout << "Failed: " << failed << '\n';

    return failed == 0 ? 0 : 1;
}
