#include "core/grid.hpp"

namespace cfd {

    Grid::Grid(std::size_t width, std::size_t height, 
        double resolution):
        width_(width), height_(height), resolution_(resolution),
        velocity_(width, height, resolution), 
        pressure_(width, height),
        bc_(width, height) {}
        
    
}