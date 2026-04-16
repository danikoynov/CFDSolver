#include "core/simulator.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <iomanip>

const double TOL = 1e-9;
const double DIV_TOL_SMALL = 1e-8;
const double DIV_TOL_BIG = 1e-8;
const int STEPS = 100;

void test_zero_state_stays_zero() {
    cfd::Simulator sim(10, 10, 0.1, 1.225, false, 0.01);

    for (int i = 0; i < STEPS; i++)
        sim.tick();

    const auto& grid = sim.grid();

    for (int i = 0; i < 10; i++)
        for (int j = 0; j < 10; j++)
            if (std::abs(grid.pressure().get_p(i, j)) > TOL)
                throw std::runtime_error(
                    "grid pressure is non-zero at cell (" +
                    std::to_string(i) + ", " + std::to_string(j) + ")"
                );

    for (int i = 0; i < 11; i++)
        for (int j = 0; j < 10; j++)
            if (std::abs(grid.velocity().get_u(i, j)) > TOL)
                throw std::runtime_error(
                    "grid u velocity component is non-zero at face (" +
                    std::to_string(i) + ", " + std::to_string(j) + ")"
                );

    for (int i = 0; i < 10; i++)
        for (int j = 0; j < 11; j++)
            if (std::abs(grid.velocity().get_v(i, j)) > TOL)
                throw std::runtime_error(
                    "grid v velocity component is non-zero at face (" +
                    std::to_string(i) + ", " + std::to_string(j) + ")"
                );
}

void test_projection_makes_small_field_divergence_free() {
    cfd::Simulator sim(3, 3, 0.1, 1.225, false, 0.01);

    auto& velocity = sim.grid().velocity();

    velocity.get_u(1, 1) = 1.0;
    velocity.get_u(2, 1) = -0.5;
    velocity.get_v(1, 1) = 0.75;
    velocity.get_v(1, 2) = -0.25;

    double max_div_before = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            max_div_before = std::max(
                max_div_before,
                std::abs(sim.grid().velocity().get_divergence(i, j))
            );
        }
    }

    if (max_div_before < 1e-12)
        throw std::runtime_error("small test setup failed: initial divergence is zero");

    sim.project(0.01);

    /*std::cout << "U Field"<<std::endl;
    for (int i = 0; i < 4; i ++, std::cout << std::endl) {
        for (int j = 0; j < 3; j ++)
            std::cout << std::setw(8) << sim.grid().velocity().get_u(i, j) << " ";
    }

    std::cout << "V Field"<<std::endl;
    for (int i = 0; i < 3; i ++, std::cout << std::endl) {
        for (int j = 0; j < 4; j ++)
            std::cout << std::setw(8) << sim.grid().velocity().get_v(i, j) << " ";
    }*/

    double max_div_after = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double div = sim.grid().velocity().get_divergence(i, j);
            max_div_after = std::max(max_div_after, std::abs(div));

            if (std::abs(div) > DIV_TOL_SMALL) {
                throw std::runtime_error(
                    "small projection test failed: divergence at cell (" +
                    std::to_string(i) + ", " +
                    std::to_string(j) + "), value = " +
                    std::to_string(div)
                );
            }
        }
    }

    if (max_div_after >= max_div_before) {
        throw std::runtime_error(
            "small projection test failed: divergence did not decrease"
        );
    }
}

void test_projection_makes_big_field_divergence_free() {
    cfd::Simulator sim(10, 10, 0.1, 1.225, false, 0.01);

    auto& velocity = sim.grid().velocity();

    velocity.get_u(5, 5) = 1.0;
    velocity.get_u(6, 5) = -1.0;
    velocity.get_v(5, 5) = 0.5;
    velocity.get_v(5, 6) = -0.5;

    double max_div_before = 0.0;
    for (int i = 0; i < static_cast<int>(sim.grid().width()); ++i) {
        for (int j = 0; j < static_cast<int>(sim.grid().height()); ++j) {
            max_div_before = std::max(
                max_div_before,
                std::abs(sim.grid().velocity().get_divergence(i, j))
            );
        }
    }

    if (max_div_before < 1e-12)
        throw std::runtime_error("big test setup failed: initial divergence is zero");

    sim.project(0.01);

    const auto& grid = sim.grid();

    double max_div_after = 0.0;
    for (int i = 0; i < static_cast<int>(grid.width()); ++i) {
        for (int j = 0; j < static_cast<int>(grid.height()); ++j) {
            double div = grid.velocity().get_divergence(i, j);
            max_div_after = std::max(max_div_after, std::abs(div));

            if (std::abs(div) > DIV_TOL_BIG) {
                throw std::runtime_error(
                    "big projection test failed: divergence at cell (" +
                    std::to_string(i) + ", " +
                    std::to_string(j) + "), value = " +
                    std::to_string(div)
                );
            }
        }
    }

    if (max_div_after >= max_div_before) {
        throw std::runtime_error(
            "big projection test failed: divergence did not decrease"
        );
    }
}

void test_body_forces_update_velocity_correctly() {
    cfd::Simulator sim(3, 3, 0.1, 1.225, false, 0.01);

    const double dt = 0.2;
    const double ax = 1.5;
    const double ay = -9.81;

    sim.apply_body_forces(dt, ax, ay);

    const auto& velocity = sim.grid().velocity();

    for (int i = 0; i <= 3; ++i)
        for (int j = 0; j < 3; ++j)
            if (std::abs(velocity.get_u(i, j) - dt * ax) > 1e-12)
                throw std::runtime_error("incorrect u update from body force");

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j <= 3; ++j)
            if (std::abs(velocity.get_v(i, j) - dt * ay) > 1e-12)
                throw std::runtime_error("incorrect v update from body force");
}
int main() {
    int passed = 0;
    int failed = 0;

    auto run = [&](void (*test)(), const char* name) {
        try {
            test();
            std::cout << "[PASS] " << name << '\n';
            passed++;
        } catch (const std::exception& e) {
            std::cout << "[FAIL] " << name << " -> " << e.what() << '\n';
            failed++;
        }
    };

    run(test_zero_state_stays_zero, "constant pressure and velocity from zero state");
    run(test_projection_makes_small_field_divergence_free, "small 3x3 divergence free projection");
    run(test_projection_makes_big_field_divergence_free, "big 10x10 divergence free projection");
    run(test_body_forces_update_velocity_correctly, "gravity applied correctly on fluid");

    std::cout << "\nPassed: " << passed << '\n';
    std::cout << "Failed: " << failed << '\n';

    return failed == 0 ? 0 : 1;
}