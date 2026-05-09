#include "core/simulator.hpp"
#include "setups/ldc.hpp"

#include <iostream>

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
    std::size_t height = 300;
    std::size_t width = 300;

    double side = 10.0;
    double resolution = side / static_cast<double>(height);
    double fluid_density = 1.225;

    double re = 100;
    double lid_velocity = 1.0;
    double viscosity = (side * lid_velocity) / re;

    bool apply_gravity = 0.0;

    cfd::Simulator sim(
        width,
        height,
        resolution,
        fluid_density,
        apply_gravity,
        viscosity
    );
    
    sim.set_profile_cgm();

    if (SETUP_TYPE == SetupType::LDC) {
        cfd::setups::LidDrivenCavity setup(lid_velocity);
        setup.impose_boundary_conditions(sim.grid());
    }

    for (unsigned int iter = 0; iter < ITERATIONS; iter ++) {
        sim.tick();
    }
}

int main() {
    profile();
    return 0;
}