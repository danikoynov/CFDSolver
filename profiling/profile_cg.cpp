#include "core/simulator.hpp"
#include "setups/ldc.hpp"

#include <iostream>
#include <filesystem>

/*
Automate running on multiple grid sizes and 
then save the data in a file to plot
*/

enum class SetupType {
    LDC
};

SetupType SETUP_TYPE = SetupType::LDC;

const unsigned int ITERATIONS = 200;

void profile() {
    namespace fs = std::filesystem;

    std::size_t height = 200;
    std::size_t width = 200;

    double side = 10.0;
    double resolution = side / static_cast<double>(height);
    double fluid_density = 1.225;

    double re = 100;
    double lid_velocity = 1.0;
    double viscosity = (side * lid_velocity) / re;

    bool apply_gravity = false;

    cfd::Simulator sim(
        width,
        height,
        resolution,
        fluid_density,
        apply_gravity,
        viscosity
    );

    fs::path source_file_path = fs::path(__FILE__);
    fs::path source_dir = source_file_path.parent_path();
    fs::path data_dir = source_dir / "data";

    sim.profile_cg_config(
        true,       // profile CG
        true,       // save CG profile data
        data_dir    // output directory
    );

    if (SETUP_TYPE == SetupType::LDC) {
        cfd::setups::LidDrivenCavity setup(lid_velocity);
        setup.impose_boundary_conditions(sim.grid());
    }

    std::cout << "Simulation setup completed" << std::endl;
    std::cout << "CG profile data directory: " << data_dir << std::endl;
    std::cout << "Begin profiling..." << std::endl;

    for (unsigned int iter = 0; iter < ITERATIONS; iter++) {
        sim.tick();
    }
}

int main() {
    profile();
    return 0;
}