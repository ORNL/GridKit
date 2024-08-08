// Created by Paul Moon 7/11/2024

#include <iostream>
#include <cmath>
#include <complex>

#include <ComponentLib/PowerFlow/Bus/BaseBus.hpp>
#include "Governor.hpp"

namespace ModelLib
{

    /*!
     * @brief Constructor for a simple governor model
     *
     * Arguments passed to ModelEvaluatorImpl:
     * - Number of equations = 6 differential + 9 algebraic = 15
     * - Quadratures and optimization parameters are not being used.
     */
    template <class ScalarT, typename IdxT>
    Governor<ScalarT, IdxT>::Governor(bus_type *bus, ScalarT P0, ScalarT Q0)
        : ModelEvaluatorImpl<ScalarT, IdxT>(18, 0, 0),
          H_(3.0),                    // Intertia constant
          D_(0.0),                    // Damping coefficient
          Xq_(0.5),                   // Machine reactance parameter
          Xd_(2.1),                   // Machine reactance parameter
          Xqp_(0.5),                  // Machine reactance parameter
          Xdp_(0.2),                  // Machine reactance parameter
          Xqpp_(0.18),                // Machine reactance parameter
          Xdpp_(0.18),                // Machine reactance parameter
          Rs_(0.0),                   // Winding resistance
          Tq0p_(0.75),                // Time constant
          Td0p_(7.0),                 // Time constant
          Tq0pp_(0.05),               // Time constant
          Td0pp_(0.035),              // Time constant
          Ef_(1.8336),                // Field winding voltage from the excitation system. An exciter would set this value
          X_(0.22),                   // Machine reactance parameter
          Xl_(0.15),                  // Machine reactance parameter
          omega_s_(0.0),              // Relative speed deviation. 0 is synchronous, -1 is 0 Hz.
          omega0_(2.0 * 60.0 * M_PI), // Nominal frequency
          S1_(0.05),                  // Saturation parameter
          S12_(0.2),                  // Saturation parameter
          R_(0.05),  // Droop constant
          T1_(0.5),  // Valve time delay
          T2_(2.5),  // Turbine numerator time constant
          T3_(7.5),  // Turbine delay
          Vmax_(1),  // Max valve position
          Vmin_(-1), // Min valve position
          Dt_(0),    // Turbine damping
          P0_(P0),  // Real power
          Q0_(Q0),  // Reactive power
          offsetGov_(15),
          bus_(bus)
    {
    }

    template <class ScalarT, typename IdxT>
    Governor<ScalarT, IdxT>::~Governor()
    {
    }

    /*!
     * @brief This function will be used to allocate sparse Jacobian matrices.
     *
     */
    template <class ScalarT, typename IdxT>
    int Governor<ScalarT, IdxT>::allocate()
    {
        tag_.resize(size_);

        return 0;
    }

    /**
     * @brief Initialization of the governor model
     *
     * Initialization equations are derived from Adam
     * Birchfield's paper on power electronics.
     *
     */
    template <class ScalarT, typename IdxT>
    int Governor<ScalarT, IdxT>::initialize()
    {
        // Compute values for saturation constants and real/imaginary voltage values
        double SB = ((1 / 0.2) * (1 / 0.2)) * ((sqrt(S1_) - sqrt(S12_)) * (sqrt(S1_) - sqrt(S12_)));
        double SA = 1 - sqrt(S1_ / SB);
        const double Vr = V() * cos(theta());
        const double Vi = V() * sin(theta());

        const double Ir = P0_ / V();
        const double Ii = 0.0;

        // Compute initial guess for the generator voltage phase
        const double delta = atan((Xq_ * Ir + Vi) / (-Xq_ * Ii + Vr));

        // Compute initial guesses for generator currents and potentials in d-q frame
        const double Id = Ir * sin(delta) - Ii * cos(delta);
        const double Iq = Ir * cos(delta) + Ii * sin(delta);

        const double Vq = Vr * cos(delta) + Vi * sin(delta) + Id * Xqpp_ + Iq * Rs_;
        const double Vd = Vr * sin(delta) - Vi * cos(delta) + Id * Rs_ - Iq * Xqpp_;

        // Compute initial guesses for flux
        const double Psiqpp = -Vd / (1 + omega_s_);
        const double Psidpp = Vq / (1 + omega_s_);

        const double Psipp = sqrt((Psiqpp * Psiqpp) + (Psidpp * Psidpp));
        const double ksat = SB * ((Psipp - SA) * (Psipp - SA));

        // Compute initial guess for Telec. Note that Pmech would normally be changed by a governor.
        const double Telec = (Psidpp - Id * Xdpp_) * Iq - (Psiqpp - Iq * Xdpp_) * Id;
        const double Pmech = Telec;

        // Compute remaining initial guesses.
        const double Psiqp = (1.0 / (1.0 + (Xqp_ - Xqpp_) / (Xqpp_ - Xl_))) * (-Psiqpp * (Xqp_ - Xl_) / (Xqpp_ - Xl_) + Iq * (Xqp_ - Xl_));

        const double Edp = Psiqp - Iq * (Xqp_ - Xl_);

        const double Psidp = -(1.0 / (1.0 + (Xdp_ - Xdpp_) / (Xdpp_ - Xl_))) * (Id * (Xdp_ - Xl_) - Psidpp * (Xdp_ - Xl_) / (Xdpp_ - Xl_));

        const double Eqp = Psidp + Id * (Xdp_ - Xl_);

        const double Pv_init = Vmin_;
        const double Ptx_init = Pmech * T3_ - T2_ * Pv_init-T3_*Dt_*omega_s_;
        const double Pref_init = P0_;

        y_[0] = delta;
        y_[1] = omega_s_;
        y_[2] = Eqp;
        y_[3] = Psidp;
        y_[4] = Psiqp;
        y_[5] = Edp;
        y_[6] = Psiqpp;
        y_[7] = Psidpp;
        y_[8] = Psipp;
        y_[9] = ksat;
        y_[10] = Vd;
        y_[11] = Vq;
        y_[12] = Telec;
        y_[13] = Id;
        y_[14] = Iq;
        y_[15] = Pmech;

        yp_[0] = 0.0;
        yp_[1] = 0.0;
        yp_[2] = 0.0;
        yp_[3] = 0.0;
        yp_[4] = 0.0;
        yp_[5] = 0.0;
        yp_[6] = 0.0;
        yp_[7] = 0.0;
        yp_[8] = 0.0;
        yp_[9] = 0.0;
        yp_[10] = 0.0;
        yp_[11] = 0.0;
        yp_[12] = 0.0;
        yp_[13] = 0.0;
        yp_[14] = 0.0;
        yp_[15] = 0.0;

        y_[offsetGov_ + 1] = Ptx_init;
        y_[offsetGov_ + 2] = Pv_init;

        yp_[offsetGov_ + 1] = 0.0;
        yp_[offsetGov_ + 2] = 0.0;

        for (int i = 0; i < size_; ++i) {
            std::cout << y_[i] << ", ";
        }
        std::cout << std::endl;

        return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <class ScalarT, typename IdxT>
    int Governor<ScalarT, IdxT>::tagDifferentiable()
    {
        tag_[0] = true;
        tag_[1] = true;
        tag_[2] = true;
        tag_[3] = true;
        tag_[4] = true;
        tag_[5] = true;
        for (IdxT i = 6; i < size_; ++i)
        {
            tag_[i] = false;
        }
        tag_[offsetGov_ + 1] = true;
        tag_[offsetGov_ + 2] = true;

        return 0;
    }

    /**
     * @brief Computes residual vector for the generator model.
     *
     * Residual equations are given by Adam Birchfield's paper
     * on power electronics.
     *
     */
    template <class ScalarT, typename IdxT>
    int Governor<ScalarT, IdxT>::evaluateResidual()
    {
        // These constants are the same as in the initialization state. Maybe move them to .hpp file?
        double SB = ((1 / 0.2) * (1 / 0.2)) * ((sqrt(S1_) - sqrt(S12_)) * (sqrt(S1_) - sqrt(S12_)));
        double SA = 1 - sqrt(S1_ / SB);
        const double Vr = V() * cos(theta());
        const double Vi = V() * sin(theta());
        double Pref = 1.0;

        // The first 6 functions are differential equations (involving yp_)
        f_[0] = -yp_[0] + omega0_ * y_[1];
        f_[1] = -yp_[1] + (1 / (2 * H_)) * ((y_[15] - D_ * y_[1]) / (1 + y_[1]) - y_[12]);
        f_[2] = -yp_[2] * Td0p_ + Ef_ - (y_[2] + (Xd_ - Xdp_) * (y_[13] + ((Xdp_ - Xdpp_) / ((Xdp_ - Xl_) * (Xdp_ - Xl_))) * (y_[2] - y_[3] - (Xdp_ - Xl_) * y_[13])) + y_[7] * y_[9]);
        f_[3] = -yp_[3] + (1 / Td0pp_) * (y_[2] - y_[3] - (Xdp_ - Xl_) * y_[13]);
        f_[4] = -yp_[4] + (1 / Tq0pp_) * (y_[5] - y_[4] + (Xqp_ - Xl_) * y_[14]);
        f_[5] = -yp_[5] + (1 / Tq0p_) * (-y_[5] + (((Xq_ - Xl_) / (Xd_ - Xl_)) * y_[6] * y_[9] + (Xq_ - Xqp_) * (y_[14] - ((Xqp_ - Xqpp_) / ((Xqp_ - Xl_) * (Xqp_ - Xl_)) * y_[5] + y_[14]) * (Xqp_ - Xl_) - y_[4])));

        // Remaining equations are all algebraic
        f_[6] = -y_[6] - y_[4] * ((Xqp_ - Xqpp_) / (Xqp_ - Xl_)) - y_[5] * ((Xqpp_ - Xl_) / (Xqp_ - Xl_));
        f_[7] = -y_[7] + y_[3] * ((Xdp_ - Xdpp_) / (Xdp_ - Xl_)) + y_[2] * ((Xdpp_ - Xl_) / (Xdp_ - Xl_));
        f_[8] = -y_[8] + sqrt((y_[7] * y_[7]) + (y_[6] * y_[6]));
        f_[9] = -y_[9] + SB * ((y_[8] - SA) * (y_[8] - SA));
        f_[10] = -y_[10] - y_[6] * (1 + y_[1]);
        f_[11] = -y_[11] + y_[7] * (1 + y_[1]);
        f_[12] = -y_[12] + y_[14] * (y_[7] - y_[13] * Xdpp_) - y_[13] * (y_[6] - y_[14] * Xdpp_);
        f_[13] = Vr * sin(y_[0]) - Vi * cos(y_[0]) + y_[13] * Rs_ - y_[14] * Xqpp_;
        f_[14] = Vr * cos(y_[0]) - Vi * sin(y_[0]) + y_[13] * Xqpp_ - y_[14] * Rs_;
        f_[15] = -y_[offsetGov_] + (1.0 / T3_) * (y_[offsetGov_ + 1] + T2_ * y_[offsetGov_ + 2]) - Dt_ * y_[1];

        P() += Pg();
        Q() += Qg();

        if (y_[offsetGov_ + 2] >= Vmax_) 
        {
            y_[offsetGov_ + 2] = Vmax_;
            if (((1.0/T1_) * (-y_[offsetGov_ + 2]) + (1.0/R_) * (Pref - y_[1])) > 0)
            {
                f_[offsetGov_+1] = -yp_[offsetGov_ + 2];
            }
            else
            {
                f_[offsetGov_+1] = -yp_[offsetGov_ + 2] + (1.0/T1_) * (-y_[offsetGov_ + 2] + (1.0/R_) * (Pref - y_[1]));
            }
        }
        else if (y_[offsetGov_ + 2] <= Vmin_) 
        {
            y_[offsetGov_ + 2] = Vmin_;
            if (((1.0/T1_) * (-y_[offsetGov_ + 2]) + (1.0/R_) * (Pref - y_[1])) > 0)
            {
                f_[offsetGov_+1] = -yp_[offsetGov_ + 2];
            }
            else
            {
                f_[offsetGov_+1] = -yp_[offsetGov_ + 2] + (1.0/T1_) * (-y_[offsetGov_ + 2] + (1.0/R_) * (Pref - y_[1]));
            }
        }
        else
        {
            f_[offsetGov_+1] = -yp_[offsetGov_ + 2] + (1.0/T1_) * (-y_[offsetGov_ + 2] + (1.0/R_) * (Pref - y_[1]));
        }
        f_[offsetGov_+2] = -yp_[offsetGov_ + 1] + y_[offsetGov_ + 2] - (1.0 / T3_) * (y_[offsetGov_ + 1] + T2_ * y_[offsetGov_ + 2]);
        //f_[offsetGov_+3] = -y_[offsetGov_] + (1.0 / T3_) * (y_[offsetGov_ + 1] + T2_ * y_[offsetGov_ + 2]) - Dt_ * y_[1];
        return 0;
    }

    /**
     * @brief Creates the Jacobian using J and Jp.
     *
     * Values were calculated by hand and verified with MATLAB.
     *
     */
    template <class ScalarT, typename IdxT>
    int Governor<ScalarT, IdxT>::evaluateJacobian()
    {
        return 0;
    }

    template <class ScalarT, typename IdxT>
    ScalarT Governor<ScalarT, IdxT>::Pg()
    {
        return y_[14]*y_[2] + y_[13]*y_[5] + (Xqp_ - Xdp_)*y_[13]*y_[14] - Rs_*(y_[13]*y_[13] + y_[14]*y_[14]);
    }

    /*y_[0] = delta;
        y_[1] = omega_s_;
        y_[2] = Eqp;
        y_[3] = Psidp;
        y_[4] = Psiqp;
        y_[5] = Edp;
        y_[6] = Psiqpp;
        y_[7] = Psidpp;
        y_[8] = Psipp;
        y_[9] = ksat;
        y_[10] = Vd;
        y_[11] = Vq;
        y_[12] = Telec;
        y_[13] = Id;
        y_[14] = Iq;
        y_[15] = Pmech;*/

/**
 * Generator reactive power Qg.
 *
 * \f[ Q_g = E_q' I_d - E_d' I_q - X_d' I_d^2 - X_q' I_q^2 \f]
 */
    template <class ScalarT, typename IdxT>
    ScalarT Governor<ScalarT, IdxT>::Qg()
    {
        return -y_[14]*y_[5] + y_[13]*y_[2] - Xdp_*y_[13]*y_[13] - Xqp_*y_[14]*y_[14];
    }

    /**
     * @brief This is not used, but is implemented to keep SystemModel.hpp happy.
     *
     */
    template <class ScalarT, typename IdxT>
    int Governor<ScalarT, IdxT>::evaluateIntegrand()
    {
        return 0;
    }

    /**
     * @brief This is not used, but is implemented to keep SystemModel.hpp happy.
     *
     */
    template <class ScalarT, typename IdxT>
    int Governor<ScalarT, IdxT>::initializeAdjoint()
    {
        return 0;
    }

    /**
     * @brief This is not used, but is implemented to keep SystemModel.hpp happy.
     *
     */
    template <class ScalarT, typename IdxT>
    int Governor<ScalarT, IdxT>::evaluateAdjointResidual()
    {
        return 0;
    }

    /**
     * @brief This is not used, but is implemented to keep SystemModel.hpp happy.
     *
     */
    template <class ScalarT, typename IdxT>
    int Governor<ScalarT, IdxT>::evaluateAdjointIntegrand()
    {
        return 0;
    }

    // Available template instantiations
    template class Governor<double, long int>;
    template class Governor<double, size_t>;

} // namespace ModelLib