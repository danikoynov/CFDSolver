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
    PressureField field(3, 2);

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
    PressureField field(4, 3);

    field.get_p(0, 0) = 1.5;
    field.get_p(2, 1) = -3.25;
    field.get_p(3, 2) = 8.0;

    require(
        field.get_p(0, 0) == 1.5,
        "test_get_p_write_and_read failed at (0, 0)"
    );
    require(
        field.get_p(2, 1) == -3.25,
        "test_get_p_write_and_read failed at (2, 1)"
    );
    require(
        field.get_p(3, 2) == 8.0,
        "test_get_p_write_and_read failed at (3, 2)"
    );
}

void test_const_get_p_reads_correct_value() {
    PressureField field(2, 2);
    field.get_p(1, 1) = 7.75;

    const PressureField& const_field = field;

    require(
        const_field.get_p(1, 1) == 7.75,
        "test_const_get_p_reads_correct_value failed"
    );
}

void test_get_p_accepts_valid_indices() {
    PressureField field(3, 4);

    field.get_p(0, 0) = 1.0;
    field.get_p(2, 3) = 2.0;
    field.get_p(1, 2) = 3.0;

    require(
        field.get_p(0, 0) == 1.0,
        "test_get_p_accepts_valid_indices failed at (0, 0)"
    );
    require(
        field.get_p(2, 3) == 2.0,
        "test_get_p_accepts_valid_indices failed at (2, 3)"
    );
    require(
        field.get_p(1, 2) == 3.0,
        "test_get_p_accepts_valid_indices failed at (1, 2)"
    );
}

void test_get_p_throws_on_negative_i() {
    PressureField field(3, 4);

    bool threw = false;

    try {
        field.get_p(-1, 0);
    } catch (const std::out_of_range&) {
        threw = true;
    }

    require(threw, "test_get_p_throws_on_negative_i failed");
}

void test_get_p_throws_on_negative_j() {
    PressureField field(3, 4);

    bool threw = false;

    try {
        field.get_p(0, -1);
    } catch (const std::out_of_range&) {
        threw = true;
    }

    require(threw, "test_get_p_throws_on_negative_j failed");
}

void test_get_p_throws_on_i_equal_width() {
    PressureField field(3, 4);

    bool threw = false;

    try {
        field.get_p(3, 0);
    } catch (const std::out_of_range&) {
        threw = true;
    }

    require(threw, "test_get_p_throws_on_i_equal_width failed");
}

void test_get_p_throws_on_j_equal_height() {
    PressureField field(3, 4);

    bool threw = false;

    try {
        field.get_p(0, 4);
    } catch (const std::out_of_range&) {
        threw = true;
    }

    require(threw, "test_get_p_throws_on_j_equal_height failed");
}

void test_const_get_p_throws_on_out_of_bounds_access() {
    PressureField field(2, 2);
    const PressureField& const_field = field;

    bool threw = false;

    try {
        const_field.get_p(2, 1);
    } catch (const std::out_of_range&) {
        threw = true;
    }

    require(threw, "test_const_get_p_throws_on_out_of_bounds_access failed");
}

void test_values_are_stored_independently() {
    PressureField field(3, 2);

    field.get_p(0, 0) = 1.0;
    field.get_p(0, 1) = 2.0;
    field.get_p(1, 0) = 3.0;
    field.get_p(1, 1) = 4.0;
    field.get_p(2, 0) = 5.0;
    field.get_p(2, 1) = 6.0;

    require(
        field.get_p(0, 0) == 1.0,
        "test_values_are_stored_independently failed at (0, 0)"
    );
    require(
        field.get_p(0, 1) == 2.0,
        "test_values_are_stored_independently failed at (0, 1)"
    );
    require(
        field.get_p(1, 0) == 3.0,
        "test_values_are_stored_independently failed at (1, 0)"
    );
    require(
        field.get_p(1, 1) == 4.0,
        "test_values_are_stored_independently failed at (1, 1)"
    );
    require(
        field.get_p(2, 0) == 5.0,
        "test_values_are_stored_independently failed at (2, 0)"
    );
    require(
        field.get_p(2, 1) == 6.0,
        "test_values_are_stored_independently failed at (2, 1)"
    );
}

void test_non_const_and_const_access_reference_same_value() {
    PressureField field(2, 3);
    field.get_p(1, 2) = 9.5;

    const PressureField& const_field = field;

    require(
        const_field.get_p(1, 2) == field.get_p(1, 2),
        "test_non_const_and_const_access_reference_same_value failed"
    );
}

}  // namespace

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

    run(test_get_p_accepts_valid_indices, "get_p accepts valid indices");
    run(test_get_p_throws_on_negative_i, "get_p throws on negative i");
    run(test_get_p_throws_on_negative_j, "get_p throws on negative j");
    run(test_get_p_throws_on_i_equal_width, "get_p throws on i equal width");
    run(test_get_p_throws_on_j_equal_height, "get_p throws on j equal height");
    run(
        test_const_get_p_throws_on_out_of_bounds_access,
        "const get_p throws on out of bounds access"
    );

    run(test_values_are_stored_independently, "values are stored independently");
    run(
        test_non_const_and_const_access_reference_same_value,
        "non-const and const access reference same value"
    );

    std::cout << "\nPassed: " << passed << '\n';
    std::cout << "Failed: " << failed << '\n';

    return failed == 0 ? 0 : 1;
}