#pragma once

#include <string>

namespace cfd::setups {
    class Airfoil {
        public:
            Airfoil(const std::string& code, double chord);

            double half_thickness(double X) const;
            double mean_camber_line(double X) const;

            const std::string& code() const;
            double chord() const;
            double m() const;
            double p() const;
            double t() const;

        private:
            static void validate_code(const std::string& code);
            static void decode(const std::string& code, double& m, double& p, double& t);

        private:
            std::string code_;
            double chord_;

            double m_; // maximum camber
            double p_; // location of maximum camber
            double t_; // thickness
    };
}