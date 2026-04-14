#include "core/boundary_conditions.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

using cfd::BoundaryConditions;
using cfd::CellType;

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
        require(
            bc.prescribed_p().empty(),
            "test_initially_no_prescribed_values: prescribed_p not empty"
        );
    }

    void test_initial_cell_types_are_fluid() {
        BoundaryConditions bc(4, 3);

        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 3; ++j) {
                require(
                    bc.type(i, j) == CellType::FLUID,
                    "test_initial_cell_types_are_fluid failed"
                );
            }
        }
    }

    void test_type_out_of_range_returns_boundary() {
        BoundaryConditions bc(4, 3);

        require(bc.type(-1, 0) == CellType::BOUNDARY, "type(-1,0) should be BOUNDARY");
        require(bc.type(4, 0) == CellType::BOUNDARY, "type(4,0) should be BOUNDARY");
        require(bc.type(0, -1) == CellType::BOUNDARY, "type(0,-1) should be BOUNDARY");
        require(bc.type(0, 3) == CellType::BOUNDARY, "type(0,3) should be BOUNDARY");
    }

    void test_set_cell_type_changes_value() {
        BoundaryConditions bc(4, 3);

        bc.set_cell_type(2, 1, CellType::SOLID);

        require(
            bc.type(2, 1) == CellType::SOLID,
            "test_set_cell_type_changes_value failed"
        );
    }

    void test_set_cell_type_out_of_range_throws() {
        BoundaryConditions bc(4, 3);

        bool threw = false;
        try {
            bc.set_cell_type(4, 0, CellType::SOLID);
        } catch (const std::out_of_range&) {
            threw = true;
        }

        require(threw, "test_set_cell_type_out_of_range_throws failed");
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

    void test_prescribe_p_value_inserts() {
        BoundaryConditions bc(4, 3);

        bc.prescribe_p_value(2, 1, 4.75);

        const auto& p = bc.prescribed_p();
        std::size_t id = static_cast<std::size_t>(2) * 3 + static_cast<std::size_t>(1);

        require(p.size() == 1, "test_prescribe_p_value_inserts: wrong map size");
        require(p.find(id) != p.end(), "test_prescribe_p_value_inserts: id not found");
        require(p.at(id) == 4.75, "test_prescribe_p_value_inserts: wrong stored value");
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

    void test_prescribe_p_duplicate_throws() {
        BoundaryConditions bc(4, 3);

        bc.prescribe_p_value(1, 2, 1.0);

        bool threw = false;
        try {
            bc.prescribe_p_value(1, 2, 2.0);
        } catch (const std::runtime_error&) {
            threw = true;
        }

        require(threw, "test_prescribe_p_duplicate_throws failed");
    }

    void test_prescribe_u_same_value_duplicate_keeps_size() {
        BoundaryConditions bc(4, 3);

        bc.prescribe_u_value(1, 1, 5.0);
        bc.prescribe_u_value(1, 1, 5.0);

        require(
            bc.prescribed_u().size() == 1,
            "test_prescribe_u_same_value_duplicate_keeps_size failed"
        );
    }

    void test_prescribe_v_same_value_duplicate_keeps_size() {
        BoundaryConditions bc(4, 3);

        bc.prescribe_v_value(1, 1, -2.0);
        bc.prescribe_v_value(1, 1, -2.0);

        require(
            bc.prescribed_v().size() == 1,
            "test_prescribe_v_same_value_duplicate_keeps_size failed"
        );
    }

    void test_prescribe_p_same_value_duplicate_keeps_size() {
        BoundaryConditions bc(4, 3);

        bc.prescribe_p_value(1, 1, 3.0);
        bc.prescribe_p_value(1, 1, 3.0);

        require(
            bc.prescribed_p().size() == 1,
            "test_prescribe_p_same_value_duplicate_keeps_size failed"
        );
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

    void test_prescribe_p_negative_i_throws() {
        BoundaryConditions bc(4, 3);

        bool threw = false;
        try {
            bc.prescribe_p_value(-1, 0, 1.0);
        } catch (const std::out_of_range&) {
            threw = true;
        }

        require(threw, "test_prescribe_p_negative_i_throws failed");
    }

    void test_prescribe_p_negative_j_throws() {
        BoundaryConditions bc(4, 3);

        bool threw = false;
        try {
            bc.prescribe_p_value(0, -1, 1.0);
        } catch (const std::out_of_range&) {
            threw = true;
        }

        require(threw, "test_prescribe_p_negative_j_throws failed");
    }

    void test_prescribe_p_i_too_large_throws() {
        BoundaryConditions bc(4, 3);

        bool threw = false;
        try {
            bc.prescribe_p_value(4, 0, 1.0);
        } catch (const std::out_of_range&) {
            threw = true;
        }

        require(threw, "test_prescribe_p_i_too_large_throws failed");
    }

    void test_prescribe_p_j_too_large_throws() {
        BoundaryConditions bc(4, 3);

        bool threw = false;
        try {
            bc.prescribe_p_value(0, 3, 1.0);
        } catch (const std::out_of_range&) {
            threw = true;
        }

        require(threw, "test_prescribe_p_j_too_large_throws failed");
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

    void test_prescribe_p_accepts_corner_cells() {
        BoundaryConditions bc(4, 3);

        bc.prescribe_p_value(0, 0, 1.0);
        bc.prescribe_p_value(3, 2, 2.0);

        require(
            bc.prescribed_p().size() == 2,
            "test_prescribe_p_accepts_corner_cells: wrong map size"
        );
    }

    void test_prescribed_u_getter_returns_value() {
        BoundaryConditions bc(4, 3);

        bc.prescribe_u_value(4, 2, 9.5);

        require(
            bc.prescribed_u(4, 2) == 9.5,
            "test_prescribed_u_getter_returns_value failed"
        );
    }

    void test_prescribed_v_getter_returns_value() {
        BoundaryConditions bc(4, 3);

        bc.prescribe_v_value(3, 3, -7.0);

        require(
            bc.prescribed_v(3, 3) == -7.0,
            "test_prescribed_v_getter_returns_value failed"
        );
    }

    void test_prescribed_p_getter_returns_value() {
        BoundaryConditions bc(4, 3);

        bc.prescribe_p_value(3, 2, 6.25);

        require(
            bc.prescribed_p(3, 2) == 6.25,
            "test_prescribed_p_getter_returns_value failed"
        );
    }

    void test_prescribed_u_missing_value_throws() {
        BoundaryConditions bc(4, 3);

        bool threw = false;
        try {
            (void)bc.prescribed_u(1, 1);
        } catch (const std::out_of_range&) {
            threw = true;
        }

        require(threw, "test_prescribed_u_missing_value_throws failed");
    }

    void test_prescribed_v_missing_value_throws() {
        BoundaryConditions bc(4, 3);

        bool threw = false;
        try {
            (void)bc.prescribed_v(1, 1);
        } catch (const std::out_of_range&) {
            threw = true;
        }

        require(threw, "test_prescribed_v_missing_value_throws failed");
    }

    void test_prescribed_p_missing_value_throws() {
        BoundaryConditions bc(4, 3);

        bool threw = false;
        try {
            (void)bc.prescribed_p(1, 1);
        } catch (const std::out_of_range&) {
            threw = true;
        }

        require(threw, "test_prescribed_p_missing_value_throws failed");
    }

    void test_is_u_prescribed_true_false() {
        BoundaryConditions bc(4, 3);

        bc.prescribe_u_value(2, 1, 3.0);

        require(bc.is_u_prescribed(2, 1), "test_is_u_prescribed_true_false: expected true");
        require(!bc.is_u_prescribed(1, 1), "test_is_u_prescribed_true_false: expected false");
    }

    void test_is_v_prescribed_true_false() {
        BoundaryConditions bc(4, 3);

        bc.prescribe_v_value(2, 1, 3.0);

        require(bc.is_v_prescribed(2, 1), "test_is_v_prescribed_true_false: expected true");
        require(!bc.is_v_prescribed(1, 1), "test_is_v_prescribed_true_false: expected false");
    }

    void test_is_p_prescribed_true_false() {
        BoundaryConditions bc(4, 3);

        bc.prescribe_p_value(2, 1, 3.0);

        require(bc.is_p_prescribed(2, 1), "test_is_p_prescribed_true_false: expected true");
        require(!bc.is_p_prescribed(1, 1), "test_is_p_prescribed_true_false: expected false");
    }

    void test_is_u_prescribed_out_of_range_throws() {
        BoundaryConditions bc(4, 3);

        bool threw = false;
        try {
            (void)bc.is_u_prescribed(6, 0);
        } catch (const std::out_of_range&) {
            threw = true;
        }

        require(threw, "test_is_u_prescribed_out_of_range_throws failed");
    }

    void test_is_v_prescribed_out_of_range_throws() {
        BoundaryConditions bc(4, 3);

        bool threw = false;
        try {
            (void)bc.is_v_prescribed(4, 0);
        } catch (const std::out_of_range&) {
            threw = true;
        }

        require(threw, "test_is_v_prescribed_out_of_range_throws failed");
    }

    void test_is_p_prescribed_out_of_range_throws() {
        BoundaryConditions bc(4, 3);

        bool threw = false;
        try {
            (void)bc.is_p_prescribed(4, 0);
        } catch (const std::out_of_range&) {
            threw = true;
        }

        require(threw, "test_is_p_prescribed_out_of_range_throws failed");
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
        require(
            bc.prescribed_p().empty(),
            "test_u_and_v_storage_are_independent: p should still be empty"
        );
    }

    void test_u_v_p_storage_are_all_independent() {
        BoundaryConditions bc(4, 3);

        bc.prescribe_u_value(1, 1, 5.0);
        bc.prescribe_v_value(1, 1, -5.0);
        bc.prescribe_p_value(1, 1, 2.5);

        require(bc.prescribed_u().size() == 1, "u size wrong");
        require(bc.prescribed_v().size() == 1, "v size wrong");
        require(bc.prescribed_p().size() == 1, "p size wrong");

        require(bc.prescribed_u(1, 1) == 5.0, "u value wrong");
        require(bc.prescribed_v(1, 1) == -5.0, "v value wrong");
        require(bc.prescribed_p(1, 1) == 2.5, "p value wrong");
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

    run(test_initial_cell_types_are_fluid, "initial cell types are fluid");
    run(test_type_out_of_range_returns_boundary, "type out of range returns boundary");
    run(test_set_cell_type_changes_value, "set cell type changes value");
    run(test_set_cell_type_out_of_range_throws, "set cell type out of range throws");

    run(test_prescribe_u_value_inserts, "prescribe u value inserts");
    run(test_prescribe_v_value_inserts, "prescribe v value inserts");
    run(test_prescribe_p_value_inserts, "prescribe p value inserts");

    run(test_prescribe_u_duplicate_throws, "prescribe u duplicate throws");
    run(test_prescribe_v_duplicate_throws, "prescribe v duplicate throws");
    run(test_prescribe_p_duplicate_throws, "prescribe p duplicate throws");

    run(test_prescribe_u_same_value_duplicate_keeps_size, "prescribe u same value duplicate keeps size");
    run(test_prescribe_v_same_value_duplicate_keeps_size, "prescribe v same value duplicate keeps size");
    run(test_prescribe_p_same_value_duplicate_keeps_size, "prescribe p same value duplicate keeps size");

    run(test_prescribe_u_negative_i_throws, "prescribe u negative i throws");
    run(test_prescribe_u_negative_j_throws, "prescribe u negative j throws");
    run(test_prescribe_u_i_too_large_throws, "prescribe u i too large throws");
    run(test_prescribe_u_j_too_large_throws, "prescribe u j too large throws");

    run(test_prescribe_v_negative_i_throws, "prescribe v negative i throws");
    run(test_prescribe_v_negative_j_throws, "prescribe v negative j throws");
    run(test_prescribe_v_i_too_large_throws, "prescribe v i too large throws");
    run(test_prescribe_v_j_too_large_throws, "prescribe v j too large throws");

    run(test_prescribe_p_negative_i_throws, "prescribe p negative i throws");
    run(test_prescribe_p_negative_j_throws, "prescribe p negative j throws");
    run(test_prescribe_p_i_too_large_throws, "prescribe p i too large throws");
    run(test_prescribe_p_j_too_large_throws, "prescribe p j too large throws");

    run(test_prescribe_u_accepts_boundary_indices, "prescribe u accepts boundary indices");
    run(test_prescribe_v_accepts_boundary_indices, "prescribe v accepts boundary indices");
    run(test_prescribe_p_accepts_corner_cells, "prescribe p accepts corner cells");

    run(test_prescribed_u_getter_returns_value, "prescribed u getter returns value");
    run(test_prescribed_v_getter_returns_value, "prescribed v getter returns value");
    run(test_prescribed_p_getter_returns_value, "prescribed p getter returns value");

    run(test_prescribed_u_missing_value_throws, "prescribed u missing value throws");
    run(test_prescribed_v_missing_value_throws, "prescribed v missing value throws");
    run(test_prescribed_p_missing_value_throws, "prescribed p missing value throws");

    run(test_is_u_prescribed_true_false, "is u prescribed true false");
    run(test_is_v_prescribed_true_false, "is v prescribed true false");
    run(test_is_p_prescribed_true_false, "is p prescribed true false");

    run(test_is_u_prescribed_out_of_range_throws, "is u prescribed out of range throws");
    run(test_is_v_prescribed_out_of_range_throws, "is v prescribed out of range throws");
    run(test_is_p_prescribed_out_of_range_throws, "is p prescribed out of range throws");

    run(test_u_and_v_storage_are_independent, "u and v storage are independent");
    run(test_u_v_p_storage_are_all_independent, "u v p storage are all independent");

    std::cout << "\nPassed: " << passed << '\n';
    std::cout << "Failed: " << failed << '\n';

    return failed == 0 ? 0 : 1;
}