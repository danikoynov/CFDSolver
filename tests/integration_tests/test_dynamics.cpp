#include "core/simulator.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

const double TOL = 1e-9;

void test_uniform_velocity_field_is_preserved_by_advection() {
    cfd::Simulator sim(5, 5, 0.1, 0.0, 1.225, false);

    auto& velocity = sim.grid().velocity();

    velocity.set_outside_velocity({2.0, -1.0});

    for (int i = 0; i <= 5; ++i)
        for (int j = 0; j < 5; ++j)
            velocity.get_u(i, j) = 2.0;

    for (int i = 0; i < 5; ++i)
        for (int j = 0; j <= 5; ++j)
            velocity.get_v(i, j) = -1.0;

    sim.advect(0.01);

    const auto& advected = sim.grid().velocity();

    for (int i = 0; i <= 5; ++i)
        for (int j = 0; j < 5; ++j)
            if (std::abs(advected.get_u(i, j) - 2.0) > TOL)
                throw std::runtime_error(
                    "u not preserved at face (" +
                    std::to_string(i) + ", " +
                    std::to_string(j) + "), value = " +
                    std::to_string(advected.get_u(i, j))
                );

    for (int i = 0; i < 5; ++i)
        for (int j = 0; j <= 5; ++j)
            if (std::abs(advected.get_v(i, j) - (-1.0)) > TOL)
                throw std::runtime_error(
                    "v not preserved at face (" +
                    std::to_string(i) + ", " +
                    std::to_string(j) + "), value = " +
                    std::to_string(advected.get_v(i, j))
                );
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

    run(
        test_uniform_velocity_field_is_preserved_by_advection,
        "uniform velocity field is preserved by advection"
    );

    std::cout << "\nPassed: " << passed << '\n';
    std::cout << "Failed: " << failed << '\n';

    return failed == 0 ? 0 : 1;
}