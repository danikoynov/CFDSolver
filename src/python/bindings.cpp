#include <pybind11/pybind11.h>
#include "core/simulator.hpp"
#include "core/grid.hpp"
#include "core/pressure_field.hpp"
#include "core/boundary_conditions.hpp"
#include "core/velocity_field.hpp"

namespace py = pybind11;

PYBIND11_MODULE(cfdsolver_py, m) {
    m.doc() = "Python bindings for CFDSolver";

    py::enum_<cfd::CellType>(m, "CellType")
        .value("FLUID", cfd::CellType::FLUID)
        .value("SOLID", cfd::CellType::SOLID)
        .value("BOUNDARY", cfd::CellType::BOUNDARY);

    py::class_<cfd::BoundaryConditions>(m, "BoundaryConditions")
        .def("prescribe_u_value", &cfd::BoundaryConditions::prescribe_u_value)
        .def("prescribe_v_value", &cfd::BoundaryConditions::prescribe_v_value)
        .def("prescribe_p_value", &cfd::BoundaryConditions::prescribe_p_value)
        .def("is_u_prescribed", &cfd::BoundaryConditions::is_u_prescribed)
        .def("is_v_prescribed", &cfd::BoundaryConditions::is_v_prescribed)
        .def("is_p_prescribed", &cfd::BoundaryConditions::is_p_prescribed)
        .def("prescribed_u", py::overload_cast<int, int>(&cfd::BoundaryConditions::prescribed_u, py::const_))
        .def("prescribed_v", py::overload_cast<int, int>(&cfd::BoundaryConditions::prescribed_v, py::const_))
        .def("prescribed_p", py::overload_cast<int, int>(&cfd::BoundaryConditions::prescribed_p, py::const_))
        .def("type", &cfd::BoundaryConditions::type)
        .def("set_cell_type", &cfd::BoundaryConditions::set_cell_type);
    
        py::class_<cfd::PressureField>(m, "PressureField")
        .def(
            "get",
            py::overload_cast<int, int>(
                &cfd::PressureField::get_p,
                py::const_
            )
        );

    py::class_<cfd::VelocityField>(m, "VelocityField")
        .def("get_u", py::overload_cast<int, int>(&cfd::VelocityField::get_u, py::const_))
        .def("get_v", py::overload_cast<int, int>(&cfd::VelocityField::get_v, py::const_));

    py::class_<cfd::Grid>(m, "Grid")
        .def("pressure",
             py::overload_cast<>(&cfd::Grid::pressure),
             py::return_value_policy::reference_internal)
        .def("boundary_conditions",
             py::overload_cast<>(&cfd::Grid::boundary_conditions),
             py::return_value_policy::reference_internal)
        .def("velocity",
             py::overload_cast<>(&cfd::Grid::velocity),
             py::return_value_policy::reference_internal)
        .def("width", &cfd::Grid::width)
        .def("height", &cfd::Grid::height);

    py::class_<cfd::Simulator>(m, "Simulator")
        .def(py::init<std::size_t, std::size_t, double, double, bool, double>())
        .def("tick", &cfd::Simulator::tick)
        .def("grid",
             py::overload_cast<>(&cfd::Simulator::grid),
             py::return_value_policy::reference_internal);
}