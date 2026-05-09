#include "setups/airfoils/airfoil.hpp"

#include <cmath>
#include <stdexcept>
#include <cctype>

namespace cfd::setups {

    Airfoil::Airfoil(const std::string& code, double chord)
        : code_(code), chord_(chord), m_(0.0), p_(0.0), t_(0.0)
    {
        validate_code(code_);

        if (chord_ <= 0.0) {
            throw std::invalid_argument("Chord must be positive");
        }

        decode(code_, m_, p_, t_);
    }

    void Airfoil::validate_code(const std::string& code)
    {
        if (code.size() != 4) {
            throw std::invalid_argument("NACA code must be 4 digits");
        }

        for (char c : code) {
            if (!std::isdigit(static_cast<unsigned char>(c))) {
                throw std::invalid_argument("NACA code must contain only digits");
            }
        }
    }

    void Airfoil::decode(const std::string& code, double& m, double& p, double& t)
    {
        m = static_cast<double>(code[0] - '0') / 100.0;
        p = static_cast<double>(code[1] - '0') / 10.0;
        t = static_cast<double>((code[2] - '0') * 10 + (code[3] - '0')) / 100.0;
    }

    double Airfoil::half_thickness(double X) const
    {
        return 5.0 * t_ * chord_ * (
            0.2969 * std::sqrt(X)
            - 0.1260 * X
            - 0.3516 * X * X
            + 0.2843 * X * X * X
            - 0.1015 * X * X * X * X
        );
    }

    double Airfoil::mean_camber_line(double X) const
    {
        // Symmetric airfoil, e.g. NACA 0012
        if (m_ == 0.0 || p_ == 0.0) {
            return 0.0;
        }

        if (X < p_) {
            return chord_ * m_ / (p_ * p_) * (
                2.0 * p_ * X - X * X
            );
        }

        return chord_ * m_ / ((1.0 - p_) * (1.0 - p_)) * (
            (1.0 - 2.0 * p_)
            + 2.0 * p_ * X
            - X * X
        );
    }

    const std::string& Airfoil::code() const
    {
        return code_;
    }

    double Airfoil::chord() const
    {
        return chord_;
    }

    double Airfoil::m() const
    {
        return m_;
    }

    double Airfoil::p() const
    {
        return p_;
    }

    double Airfoil::t() const
    {
        return t_;
    }

}