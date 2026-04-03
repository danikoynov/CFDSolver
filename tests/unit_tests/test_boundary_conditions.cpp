#include "core/boundary_conditions.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

using cfd::BoundaryConditions;

namespace {

    void require(bool condition, const std::string& message) {
        if (!condition) {
            throw std::runtime_error(message);
        }
    }

    void test_initially_no_prescribed_values() {
        BoundaryConditions bc(4, 3);

        require(
            bc.prescribed_u().empty(),
            "test_initially_no_prescribed_values: prescribed_u not empty"
        );
        require(
            bc.prescribed_v().empty(),
            "test_initially_no_prescribed_values: prescribed_v not empty"
        );
    }

    void test_prescribe_u_value_inserts() {
        BoundaryConditions bc(4, 3);

        bc.prescribe_u_value(2, 1, 7.5);

        const auto& u = bc.prescribed_u();
        std::size_t id = static_cast<std::size_t>(2) * 3 + static_cast<std::size_t>(1);

        require(u.size() == 1, "test_prescribe_u_value_inserts: wrong map size");
        require(u.find(id) != u.end(), "test_prescribe_u_value_inserts: id not found");
        require(u.at(id) == 7.5, "test_prescribe_u_value_inserts: wrong stored value");
    }

    void test_prescribe_v_value_inserts() {
        BoundaryConditions bc(4, 3);

        bc.prescribe_v_value(2, 1, -3.25);

        const auto& v = bc.prescribed_v();
        std::size_t id = static_cast<std::size_t>(2) * (3 + 1) + static_cast<std::size_t>(1);

        require(v.size() == 1, "test_prescribe_v_value_inserts: wrong map size");
        require(v.find(id) != v.end(), "test_prescribe_v_value_inserts: id not found");
        require(v.at(id) == -3.25, "test_prescribe_v_value_inserts: wrong stored value");
    }

    void test_prescribe_u_duplicate_throws() {
        BoundaryConditions bc(4, 3);

        bc.prescribe_u_value(1, 2, 1.0);

        bool threw = false;
        try {
            bc.prescribe_u_value(1, 2, 2.0);
        } catch (const std::runtime_error&) {
            threw = true;
        }

        require(threw, "test_prescribe_u_duplicate_throws failed");
    }

    void test_prescribe_v_duplicate_throws() {
        BoundaryConditions bc(4, 3);

        bc.prescribe_v_value(1, 2, 1.0);

        bool threw = false;
        try {
            bc.prescribe_v_value(1, 2, 2.0);
        } catch (const std::runtime_error&) {
            threw = true;
        }

        require(threw, "test_prescribe_v_duplicate_throws failed");
    }

    void test_prescribe_u_negative_i_throws() {
        BoundaryConditions bc(4, 3);

        bool threw = false;
        try {
            bc.prescribe_u_value(-1, 0, 1.0);
        } catch (const std::out_of_range&) {
            threw = true;
        }

        require(threw, "test_prescribe_u_negative_i_throws failed");
    }

    void test_prescribe_u_negative_j_throws() {
        BoundaryConditions bc(4, 3);

        bool threw = false;
        try {
            bc.prescribe_u_value(0, -1, 1.0);
        } catch (const std::out_of_range&) {
            threw = true;
        }

        require(threw, "test_prescribe_u_negative_j_throws failed");
    }

    void test_prescribe_u_i_too_large_throws() {
        BoundaryConditions bc(4, 3);

        bool threw = false;
        try {
            bc.prescribe_u_value(5, 0, 1.0);
        } catch (const std::out_of_range&) {
            threw = true;
        }

        require(threw, "test_prescribe_u_i_too_large_throws failed");
    }

    void test_prescribe_u_j_too_large_throws() {
        BoundaryConditions bc(4, 3);

        bool threw = false;
        try {
            bc.prescribe_u_value(0, 3, 1.0);
        } catch (const std::out_of_range&) {
            threw = true;
        }

        require(threw, "test_prescribe_u_j_too_large_throws failed");
    }

    void test_prescribe_v_negative_i_throws() {
        BoundaryConditions bc(4, 3);

        bool threw = false;
        try {
            bc.prescribe_v_value(-1, 0, 1.0);
        } catch (const std::out_of_range&) {
            threw = true;
        }

        require(threw, "test_prescribe_v_negative_i_throws failed");
    }

    void test_prescribe_v_negative_j_throws() {
        BoundaryConditions bc(4, 3);

        bool threw = false;
        try {
            bc.prescribe_v_value(0, -1, 1.0);
        } catch (const std::out_of_range&) {
            threw = true;
        }

        require(threw, "test_prescribe_v_negative_j_throws failed");
    }

    void test_prescribe_v_i_too_large_throws() {
        BoundaryConditions bc(4, 3);

        bool threw = false;
        try {
            bc.prescribe_v_value(4, 0, 1.0);
        } catch (const std::out_of_range&) {
            threw = true;
        }

        require(threw, "test_prescribe_v_i_too_large_throws failed");
    }

    void test_prescribe_v_j_too_large_throws() {
        BoundaryConditions bc(4, 3);

        bool threw = false;
        try {
            bc.prescribe_v_value(0, 4, 1.0);
        } catch (const std::out_of_range&) {
            threw = true;
        }

        require(threw, "test_prescribe_v_j_too_large_throws failed");
    }

    void test_prescribe_u_accepts_boundary_indices() {
        BoundaryConditions bc(4, 3);

        bc.prescribe_u_value(0, 0, 1.0);
        bc.prescribe_u_value(4, 2, 2.0);

        require(
            bc.prescribed_u().size() == 2,
            "test_prescribe_u_accepts_boundary_indices: wrong map size"
        );
    }

    void test_prescribe_v_accepts_boundary_indices() {
        BoundaryConditions bc(4, 3);

        bc.prescribe_v_value(0, 0, 1.0);
        bc.prescribe_v_value(3, 3, 2.0);

        require(
            bc.prescribed_v().size() == 2,
            "test_prescribe_v_accepts_boundary_indices: wrong map size"
        );
    }

    void test_u_and_v_storage_are_independent() {
        BoundaryConditions bc(4, 3);

        bc.prescribe_u_value(1, 1, 5.0);
        bc.prescribe_v_value(1, 1, -5.0);

        require(
            bc.prescribed_u().size() == 1,
            "test_u_and_v_storage_are_independent: wrong u size"
        );
        require(
            bc.prescribed_v().size() == 1,
            "test_u_and_v_storage_are_independent: wrong v size"
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

    run(test_initially_no_prescribed_values, "initially no prescribed values");

    run(test_prescribe_u_value_inserts, "prescribe u value inserts");
    run(test_prescribe_v_value_inserts, "prescribe v value inserts");

    run(test_prescribe_u_duplicate_throws, "prescribe u duplicate throws");
    run(test_prescribe_v_duplicate_throws, "prescribe v duplicate throws");

    run(test_prescribe_u_negative_i_throws, "prescribe u negative i throws");
    run(test_prescribe_u_negative_j_throws, "prescribe u negative j throws");
    run(test_prescribe_u_i_too_large_throws, "prescribe u i too large throws");
    run(test_prescribe_u_j_too_large_throws, "prescribe u j too large throws");

    run(test_prescribe_v_negative_i_throws, "prescribe v negative i throws");
    run(test_prescribe_v_negative_j_throws, "prescribe v negative j throws");
    run(test_prescribe_v_i_too_large_throws, "prescribe v i too large throws");
    run(test_prescribe_v_j_too_large_throws, "prescribe v j too large throws");

    run(test_prescribe_u_accepts_boundary_indices, "prescribe u accepts boundary indices");
    run(test_prescribe_v_accepts_boundary_indices, "prescribe v accepts boundary indices");

    run(test_u_and_v_storage_are_independent, "u and v storage are independent");

    std::cout << "\nPassed: " << passed << '\n';
    std::cout << "Failed: " << failed << '\n';

    return failed == 0 ? 0 : 1;
}