#include "core/velocity_field.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using cfd::Velocity2D;
using cfd::VelocityField;

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

void require_velocity_close(
    const Velocity2D& v,
    double expected_u,
    double expected_v,
    const std::string& test_name
) {
    if (!approx_equal(v.u, expected_u)) {
        throw std::runtime_error(test_name + ": wrong u component");
    }
    if (!approx_equal(v.v, expected_v)) {
        throw std::runtime_error(test_name + ": wrong v component");
    }
}

void test_get_u_set_and_get() {
    VelocityField field(3, 2, 1.0);

    field.get_u(0, 0) = 1.5;
    field.get_u(3, 1) = -2.0;

    require(
        approx_equal(field.get_u(0, 0), 1.5),
        "test_get_u_set_and_get failed at (0,0)"
    );
    require(
        approx_equal(field.get_u(3, 1), -2.0),
        "test_get_u_set_and_get failed at (3,1)"
    );
}

void test_get_v_set_and_get() {
    VelocityField field(3, 2, 1.0);

    field.get_v(0, 0) = 2.5;
    field.get_v(2, 2) = -1.25;

    require(
        approx_equal(field.get_v(0, 0), 2.5),
        "test_get_v_set_and_get failed at (0,0)"
    );
    require(
        approx_equal(field.get_v(2, 2), -1.25),
        "test_get_v_set_and_get failed at (2,2)"
    );
}

void test_get_u_out_of_bounds_throws() {
    VelocityField field(3, 2, 1.0);

    bool threw = false;
    try {
        field.get_u(4, 0);
    } catch (const std::out_of_range&) {
        threw = true;
    }

    require(threw, "test_get_u_out_of_bounds_throws failed");
}

void test_get_v_out_of_bounds_throws() {
    VelocityField field(3, 2, 1.0);

    bool threw = false;
    try {
        field.get_v(3, 0);
    } catch (const std::out_of_range&) {
        threw = true;
    }

    require(threw, "test_get_v_out_of_bounds_throws failed");
}

void test_sample_at_vertical_face_uniform_field() {
    VelocityField field(3, 3, 1.0);

    for (int i = 0; i <= 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            field.get_u(i, j) = 2.0;
        }
    }

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j <= 3; ++j) {
            field.get_v(i, j) = -1.0;
        }
    }

    Velocity2D vel = field.sample_at_vertical_face(1, 1);

    require_velocity_close(
        vel,
        2.0,
        -1.0,
        "test_sample_at_vertical_face_uniform_field"
    );
}

void test_sample_at_horizontal_face_uniform_field() {
    VelocityField field(3, 3, 1.0);

    for (int i = 0; i <= 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            field.get_u(i, j) = 2.0;
        }
    }

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j <= 3; ++j) {
            field.get_v(i, j) = -1.0;
        }
    }

    Velocity2D vel = field.sample_at_horizontal_face(1, 1);

    require_velocity_close(
        vel,
        2.0,
        -1.0,
        "test_sample_at_horizontal_face_uniform_field"
    );
}

void test_sample_at_corner_uniform_field() {
    VelocityField field(3, 3, 1.0);

    for (int i = 0; i <= 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            field.get_u(i, j) = 4.0;
        }
    }

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j <= 3; ++j) {
            field.get_v(i, j) = 7.0;
        }
    }

    Velocity2D vel = field.sample_at_corner(1, 1);

    require_velocity_close(
        vel,
        4.0,
        7.0,
        "test_sample_at_corner_uniform_field"
    );
}

void test_sample_at_coordinates_uniform_field() {
    VelocityField field(4, 4, 1.0);

    for (int i = 0; i <= 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            field.get_u(i, j) = 3.0;
        }
    }

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j <= 4; ++j) {
            field.get_v(i, j) = -2.0;
        }
    }

    Velocity2D vel = field.sample_at_coordinates(1.3, 1.4);

    require_velocity_close(
        vel,
        3.0,
        -2.0,
        "test_sample_at_coordinates_uniform_field"
    );
}

void test_sample_at_coordinates_clamps_to_boundary() {
    VelocityField field(3, 3, 1.0);

    for (int i = 0; i <= 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            field.get_u(i, j) = 5.0;
        }
    }

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j <= 3; ++j) {
            field.get_v(i, j) = -2.0;
        }
    }

    Velocity2D vel = field.sample_at_coordinates(-0.1, 1.0);

    require_velocity_close(
        vel,
        5.0,
        -2.0,
        "test_sample_at_coordinates_clamps_to_boundary"
    );
}

void test_sample_at_coordinates_exact_grid_corner() {
    VelocityField field(2, 2, 1.0);

    for (int i = 0; i <= 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            field.get_u(i, j) = 5.0;
        }
    }

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j <= 2; ++j) {
            field.get_v(i, j) = -3.0;
        }
    }

    Velocity2D vel = field.sample_at_coordinates(1.0, 1.0);

    require_velocity_close(
        vel,
        5.0,
        -3.0,
        "test_sample_at_coordinates_exact_grid_corner"
    );
}

void test_linear_field_reproduced_inside_cell() {
    VelocityField field(3, 3, 1.0);

    for (int i = 0; i <= 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            field.get_u(i, j) = static_cast<double>(i);
        }
    }

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j <= 3; ++j) {
            field.get_v(i, j) = static_cast<double>(j);
        }
    }

    Velocity2D vel = field.sample_at_coordinates(1.25, 1.75);

    require(
        approx_equal(vel.u, 1.25),
        "test_linear_field_reproduced_inside_cell failed for u"
    );
    require(
        approx_equal(vel.v, 1.75),
        "test_linear_field_reproduced_inside_cell failed for v"
    );
}

void test_boundary_face_clamping_behavior() {
    VelocityField field(2, 2, 1.0);

    for (int i = 0; i <= 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            field.get_u(i, j) = 1.0;
        }
    }

    field.get_v(0, 0) = 10.0;
    field.get_v(0, 1) = 20.0;
    field.get_v(1, 0) = 30.0;
    field.get_v(1, 1) = 40.0;
    field.get_v(0, 2) = 50.0;
    field.get_v(1, 2) = 60.0;

    Velocity2D vel = field.sample_at_vertical_face(0, 0);

    require(
        approx_equal(vel.u, 1.0),
        "test_boundary_face_clamping_behavior failed for u"
    );
    require(
        approx_equal(vel.v, 15.0),
        "test_boundary_face_clamping_behavior failed for v"
    );
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

    run(test_get_u_set_and_get, "get_u set and get");
    run(test_get_v_set_and_get, "get_v set and get");
    run(test_get_u_out_of_bounds_throws, "get_u out of bounds throws");
    run(test_get_v_out_of_bounds_throws, "get_v out of bounds throws");

    run(test_sample_at_vertical_face_uniform_field, "sample at vertical face uniform field");
    run(test_sample_at_horizontal_face_uniform_field, "sample at horizontal face uniform field");
    run(test_sample_at_corner_uniform_field, "sample at corner uniform field");

    run(test_sample_at_coordinates_uniform_field, "sample at coordinates uniform field");
    run(test_sample_at_coordinates_clamps_to_boundary, "sample outside clamps to boundary");
    run(test_sample_at_coordinates_exact_grid_corner, "sample at coordinates exact grid corner");

    run(test_linear_field_reproduced_inside_cell, "linear field reproduced inside cell");
    run(test_boundary_face_clamping_behavior, "boundary face clamping behavior");

    std::cout << "\nPassed: " << passed << '\n';
    std::cout << "Failed: " << failed << '\n';

    return failed == 0 ? 0 : 1;
}