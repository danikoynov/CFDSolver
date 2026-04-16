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
        
        const double DEFAULT_PRESSURE = 0.0;
        
        private:
            std::size_t width_, height_;
            std::vector<CellType> type_; /// Default assigned to fluid
            std::unordered_map<std::size_t, double> prescribed_u_;
            std::unordered_map<std::size_t, double> prescribed_v_;
            std::unordered_map<int, double> prescribed_p_;

        public:
            BoundaryConditions(std::size_t width, std::size_t height);
    
            void set_cell_type(int i, int j, CellType cell_type); 
            void prescribe_u_value(int i, int j, double val);
            void prescribe_v_value(int i, int j, double val);
            void prescribe_p_value(int i, int j, double val);

            bool is_u_prescribed(int i, int j) const;
            bool is_v_prescribed(int i, int j) const;
            bool is_p_prescribed(int i, int j) const;

            const std::unordered_map<std::size_t, double>& prescribed_u() const;
            const std::unordered_map<std::size_t, double>& prescribed_v() const;
            const std::unordered_map<int, double>& prescribed_p() const;

            double prescribed_u(int i, int j) const;
            double prescribed_v(int i, int j) const;
            double prescribed_p(int i, int j) const;

            CellType type(int i, int j) const;

            bool is_solid(int i, int j) const;
    };
}