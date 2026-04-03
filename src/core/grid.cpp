#include "core/grid.hpp"

namespace cfd {

    Grid::Grid(std::size_t width, std::size_t height, 
        double resolution, double outside_pressure):
        width_(width), height_(height),
        velocity_(width, height, resolution), 
        pressure_(width, height, outside_pressure),
        bc_(width, height) {}
        
    
}