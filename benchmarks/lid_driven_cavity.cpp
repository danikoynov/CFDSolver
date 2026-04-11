#include "core/simulator.hpp"
#include "core/velocity_field.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include <filesystem>
#include <unordered_map>
#include <cmath>
#include <sstream>
#include <iomanip>

constexpr bool DEBUG_PRINT = false;

namespace fs = std::filesystem;

const double CONVERGENCE_TOL = 1e-3;

struct Entry {
    double y, val;

};

std::unordered_map<int, std::vector<Entry>> read_table(
    const std::string& filename) {

    std::string subfolder = "benchmarks/reference";
    fs::path path = subfolder + "/" + filename;

    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("Could not open: " + filename);
    }
    else {
        std::cout<<"Opened: "<<filename<<std::endl;
    }

    std::unordered_map<int, std::vector<Entry>> entries;
    std::unordered_map<int, int> col_to_re;

    std::string line;
    bool first_line = true;

    
    while(std::getline(file, line)) {
        std::stringstream ss(line);
        std::string cell;
        std::vector<std::string> row;

        while(std::getline(ss, cell, ',')) {
            row.push_back(cell);
        }

        if (first_line) {
            std::size_t cols = row.size();
            for (int i = 2; i < cols; i ++) {
                std::size_t re_length = row[i].size();

                std::string re_str = row[i].substr(2, re_length - 2);
                int re_num = static_cast<int>(stod(re_str));
                
                col_to_re[i] = re_num;
            }    

            first_line = false;
        }
        else {
            std::size_t cols = row.size();
            double y = std::stod(row[1]);

            for (int i = 2; i < cols; i ++) {
                double val = std::stod(row[i]);
                int re = col_to_re[i];
                entries[re].push_back({y, val});
            }
        }
    }
    file.close();

    return entries;
}

bool has_converged(std::size_t w, std::size_t h,
    const cfd::VelocityField& A, const cfd::VelocityField& B) {
    
    for (int i = 0; i <= w; i ++)
        for (int j = 0; j < h; j ++) {
            if (std::abs(A.get_u(i, j) - B.get_u(i, j)) > CONVERGENCE_TOL) {
                if (DEBUG_PRINT) {
                    std::cout<<"Diff: "<<A.get_u(i, j) - B.get_u(i, j)<< 
                        " at: "<<i<<" "<<j<<" horizontal"<<std::endl;
                }
                return false;
            }
        }

    for (int i = 0; i < w; i ++)
        for (int j = 0; j <= h; j ++) {
            if (std::abs(A.get_v(i, j) - B.get_v(i, j)) > CONVERGENCE_TOL) {
                if (DEBUG_PRINT) {
                    std::cout<<"Diff: "<<A.get_v(i, j) - B.get_v(i, j)<< 
                        " at: "<<i<<" "<<j<<" vertical"<<std::endl;   
                } 
                return false;
            }
        }

    return true;
}

void set_to_solid(int i, int j, cfd::BoundaryConditions& bc) {
    bc.set_cell_type(i, j, cfd::CellType::SOLID);

    bc.prescribe_u_value(i, j, 0);
    bc.prescribe_u_value(i + 1, j, 0);

    bc.prescribe_v_value(i, j, 0);
    bc.prescribe_v_value(i, j + 1, 0);
}

void apply_lid_driven_bc(cfd::Simulator& sim, double lid_velocity) {
    int w = sim.grid().width();
    int h = sim.grid().height();

    auto& bc = sim.grid().boundary_conditions();
    for (int i = 0; i < w; i ++) {
        set_to_solid(i, 0, bc);
        set_to_solid(i, h - 1, bc);
    }
    
    for (int j = 0; j < h; j ++) {
        set_to_solid(0, j, bc);
        set_to_solid(w - 1, j, bc);
    }

    for (int i = 2; i < w - 1; i ++) {
        bc.prescribe_u_value(i, h - 2, lid_velocity);
        bc.prescribe_u_value(i, 1, 0);
    }

    for (int j = 2; j < h - 1; j ++) {
        bc.prescribe_v_value(1, j, 0);
        bc.prescribe_v_value(w - 2, j, 0);
    }
}

double sample_u_vertical_centerline(const cfd::VelocityField& vf,
    double y, std::size_t w, std::size_t h, double dx) {

    if (w % 2 != 0) {
        throw std::runtime_error("Vertical centerline sampling requires even width");
    }

    int i = static_cast<int>(w / 2);

    double y_min = 0.5 * dx;
    double y_max = (static_cast<double>(h) - 0.5) * dx;

    if (y <= y_min) return vf.get_u(i, 0);
    if (y >= y_max) return vf.get_u(i, static_cast<int>(h) - 1);

    double pos = y / dx - 0.5;
    int j0 = static_cast<int>(std::floor(pos));
    int j1 = j0 + 1;

    double y0 = (static_cast<double>(j0) + 0.5) * dx;
    double y1 = (static_cast<double>(j1) + 0.5) * dx;

    double u0 = vf.get_u(i, j0);
    double u1 = vf.get_u(i, j1);

    double t = (y - y0) / (y1 - y0);
    return (1.0 - t) * u0 + t * u1;
}

double evaluate_u_error(const cfd::VelocityField& vf,
    const std::vector<Entry>& entries, 
    std::size_t h, std::size_t w, double dx) {

    double sum = 0.0;
    int n = 0;

    for (const auto& [y, val] : entries) {
        double x_cor = 0.5;
        double y_cor = dx + dx / 2 + y * (h - 3) * dx;

        double vc = sample_u_vertical_centerline(vf, y_cor, w, h, dx);

        double diff = (val - vc);
        sum = sum + diff * diff;
        n ++;
        if (DEBUG_PRINT) {
            std::cout<<"Horizontal | at "<<x_cor<<" "<<y_cor<< ": " << vc << " // " << val <<std::endl;
        }
    }

    if (DEBUG_PRINT) {
        std::cout<<"U Velocity Field"<<std::endl;
        for (int j = h - 1; j >= 0; j --) {
            for (int i = 0; i <= w; i ++) {
                std::cout<<std::setw(6)<<std::fixed<<std::setprecision(3)<<vf.get_u(i, j) << " ";
            }
            std::cout << std::endl;
        }
    }

    double rmse = std::sqrt(sum / static_cast<double>(n));
    return rmse;
}
    

double sample_v_horizontal_centerline(const cfd::VelocityField& vf,
    double x, std::size_t w, std::size_t h, double dx) {

    if (h % 2 != 0) {
        throw std::runtime_error("Horizontal centerline sampling requires even height");
    }

    int j = static_cast<int>(h / 2);

    double x_min = 0.5 * dx;
    double x_max = (static_cast<double>(w) - 0.5) * dx;

    if (x <= x_min) return vf.get_v(0, j);
    if (x >= x_max) return vf.get_v(static_cast<int>(w) - 1, j);

    double pos = x / dx - 0.5;
    int i0 = static_cast<int>(std::floor(pos));
    int i1 = i0 + 1;

    double x0 = (static_cast<double>(i0) + 0.5) * dx;
    double x1 = (static_cast<double>(i1) + 0.5) * dx;

    double v0 = vf.get_v(i0, j);
    double v1 = vf.get_v(i1, j);

    double t = (x - x0) / (x1 - x0);
    return (1.0 - t) * v0 + t * v1;
}

double evaluate_v_error(const cfd::VelocityField& vf,
    const std::vector<Entry>& entries,
    std::size_t h, std::size_t w, double dx) {

    double sum = 0.0;
    int n = 0;

    for (const auto& [x, val] : entries) {
        double x_cor = dx + dx / 2 + x * (w - 3) * dx;
        double y_cor = 0.5;

        double vc = sample_v_horizontal_centerline(vf, x_cor, w, h, dx);

        double diff = (val - vc);
        sum += diff * diff;
        n++;

        if (DEBUG_PRINT) { 
            std::cout << "Vertical | at " << x_cor << " " << y_cor
                  << ": " << vc << " // " << val << std::endl;
        }
    }
    
    if (DEBUG_PRINT) { 
        std::cout<<"V Velocity Field"<<std::endl;
        for (int j = h; j >= 0; j --) {
            for (int i = 0; i < w; i ++) {
                std::cout<<std::setw(6)<<std::fixed<<std::setprecision(3)<<vf.get_v(i, j) << " ";
            }
            std::cout << std::endl;
        }
    }

    double rmse = std::sqrt(sum / static_cast<double>(n));
    return rmse;
}


void run_test(double re, 
    const std::vector<Entry>& u_entries, const std::vector<Entry>& v_entries) {
    
    std::size_t width = 30, height = 30;
            
    double len = 1.0; /// Lenght of side
    double vel = 1.0; /// Lid velocity
            
    double viscosity = (len * vel) / re;
    double density = 1.225;
    double dx = len / static_cast<double>(height);

    cfd::Simulator sim = cfd::Simulator(width, height,
        dx, 0, density, false, viscosity);
    
    apply_lid_driven_bc(sim, vel);

    auto last_vf = sim.grid().velocity();

    for (int i = 0; true; i ++) {
        sim.tick();
        const auto& new_vf = sim.grid().velocity();
        if (has_converged(width, height, last_vf, new_vf))
            break;

        last_vf = new_vf;

        if (i % 10 == 0) {
            std::cout<<"Tick: "<<i<<std::endl;
        }

    }

    double u_rmse = evaluate_u_error(last_vf, u_entries, height, width, dx); // fix the hor or ver
    double v_rmse = evaluate_v_error(last_vf, v_entries, height, width, dx); // fix the hor or ver
    
    std::cout<<"Benchmark finished"<<std::endl;
    std::cout<<"U RMSE: "<<u_rmse<<std::endl;
    std::cout<<"V RMSE: "<<v_rmse<<std::endl;
}


int main() {
    std::string u_file = "ghia_table1_u_vertical_centerline.csv";
    std::string v_file = "ghia_table1_v_horizontal_centerline.csv";
    
    auto u_entries = read_table(u_file);
    auto v_entries = read_table(v_file);
    
    std::vector<int> re_list;
    for (const auto& [key, value]: u_entries) {
        re_list.push_back(key);
    }

    int re = 100;
    run_test(static_cast<double>(re), u_entries[re], v_entries[re]);


}