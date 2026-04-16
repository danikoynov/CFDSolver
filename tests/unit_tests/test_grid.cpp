#include "core/grid.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

using cfd::Grid;

namespace {

    void require(bool condition, const std::string& message) {
        if (!condition) {
            throw std::runtime_error(message);
        }
    }

    void test_grid_dimensions() {
        Grid grid(5, 4, 1.0);

        require(
            grid.width() == 5,
            "test_grid_dimensions failed for width"
        );
        require(
            grid.height() == 4,
            "test_grid_dimensions failed for height"
        );
    }

    void test_velocity_accessor_is_usable() {
        Grid grid(3, 2, 1.0);

        grid.velocity().get_u(0, 0) = 2.5;
        grid.velocity().get_v(1, 1) = -1.0;

        require(
            grid.velocity().get_u(0, 0) == 2.5,
            "test_velocity_accessor_is_usable failed for u"
        );
        require(
            grid.velocity().get_v(1, 1) == -1.0,
            "test_velocity_accessor_is_usable failed for v"
        );
    }

    void test_pressure_accessor_is_usable() {
        Grid grid(4, 3, 1.0);

        grid.pressure().get_p(1, 2) = 42.0;

        require(
            grid.pressure().get_p(1, 2) == 42.0,
            "test_pressure_accessor_is_usable failed"
        );
    }

    void test_boundary_conditions_accessor_is_usable() {
        Grid grid(3, 3, 1.0);

        grid.boundary_conditions().prescribe_u_value(0, 1, 5.0);

        require(
            true,
            "test_boundary_conditions_accessor_is_usable failed"
        );
    }

    void test_const_accessors_are_usable() {
        Grid grid(2, 2, 1.0);

        grid.velocity().get_u(0, 0) = 3.0;
        grid.pressure().get_p(1, 1) = 7.0;

        const Grid& const_grid = grid;

        require(
            const_grid.velocity().get_u(0, 0) == 3.0,
            "test_const_accessors_are_usable failed for velocity"
        );
        require(
            const_grid.pressure().get_p(1, 1) == 7.0,
            "test_const_accessors_are_usable failed for pressure"
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

    run(test_grid_dimensions, "grid dimensions");
    run(test_velocity_accessor_is_usable, "velocity accessor is usable");
    run(test_pressure_accessor_is_usable, "pressure accessor is usable");
    run(test_boundary_conditions_accessor_is_usable, "boundary conditions accessor is usable");
    run(test_const_accessors_are_usable, "const accessors are usable");

    std::cout << "\nPassed: " << passed << '\n';
    std::cout << "Failed: " << failed << '\n';

    return failed == 0 ? 0 : 1;
}
