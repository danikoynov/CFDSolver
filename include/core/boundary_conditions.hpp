#pragma once
#include <vector>
#include <unordered_map>
namespace cfd {

    enum class CellType{
        FLUID,
        SOLID,
        BOUNDARY
    };

    class BoundaryConditions {
        /// Represents the predescribed boundary conditions in the simulation
        
        private:
            std::size_t width_, height_;
            std::vector<CellType> type_; /// Default assigned to fluid
            std::unordered_map<std::size_t, double> prescribed_u_;
            std::unordered_map<std::size_t, double> prescribed_v_;

        public:
            BoundaryConditions(std::size_t width, std::size_t height);

            void prescribe_u_value(int i, int j, double val);
            void prescribe_v_value(int i, int j, double val);

            const std::unordered_map<std::size_t, double>& prescribed_u() const;
            const std::unordered_map<std::size_t, double>& prescribed_v() const;
    };
}