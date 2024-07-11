// Made by Paul Moon 6/7/2024

#include <iostream>
#include <cmath>
#include <complex>

#include <ComponentLib/PowerFlow/Bus/BaseBus.hpp>
#include "GENROU.hpp"

namespace ModelLib
{

    /*!
     * @brief Constructor for a simple generator model
     *
     * Arguments passed to ModelEvaluatorImpl:
     * - Number of equations = 6 differential + 9 algebraic = 15
     * - Quadratures and optimization parameters are not being used.
     */
    template <class ScalarT, typename IdxT>
    GENROU<ScalarT, IdxT>::GENROU(bus_type *bus, ScalarT P0, ScalarT Q0)
        : ModelEvaluatorImpl<ScalarT, IdxT>(15, 0, 0),
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
          Pm_(1.0),                   // Mechanical power from the prime mover. A governor would set this value
          X_(0.22),                   // Machine reactance parameter
          Xl_(0.15),                  // Machine reactance parameter
          omega_s_(0.0),              // Relative speed deviation. 0 is synchronous, -1 is 0 Hz.
          omega0_(2.0 * 60.0 * M_PI), // Nominal frequency
          S1_(0.05),                  // Saturation parameter
          S12_(0.2),                  // Saturation parameter
          P0_(1.0),                   // Real power
          Q0_(0.0),                   // Reactive power
          bus_(bus)
    {
    }

    template <class ScalarT, typename IdxT>
    GENROU<ScalarT, IdxT>::~GENROU()
    {
    }

    /*!
     * @brief This function will be used to allocate sparse Jacobian matrices.
     *
     */
    template <class ScalarT, typename IdxT>
    int GENROU<ScalarT, IdxT>::allocate()
    {
        tag_.resize(size_);

        return 0;
    }

    /**
     * @brief Initialization of the generator model
     *
     * Initialization equations are derived from Adam
     * Birchfield's paper on power electronics.
     *
     */
    template <class ScalarT, typename IdxT>
    int GENROU<ScalarT, IdxT>::initialize()
    {
        // Compute values for saturation constants and real/imaginary voltage values
        double SB = ((1 / 0.2) * (1 / 0.2)) * ((sqrt(S1_) - sqrt(S12_)) * (sqrt(S1_) - sqrt(S12_)));
        double SA = 1 - sqrt(S1_ / SB);
        const double Vr = V() * cos(theta());
        const double Vi = V() * sin(theta());

        const double Ir = P0_ / V();
        const double Ii = 0;

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

        return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <class ScalarT, typename IdxT>
    int GENROU<ScalarT, IdxT>::tagDifferentiable()
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
    int GENROU<ScalarT, IdxT>::evaluateResidual()
    {
        // These constants are the same as in the initialization state. Maybe move them to .hpp file?
        double SB = ((1 / 0.2) * (1 / 0.2)) * ((sqrt(S1_) - sqrt(S12_)) * (sqrt(S1_) - sqrt(S12_)));
        double SA = 1 - sqrt(S1_ / SB);
        const double Vr = V() * cos(theta());
        const double Vi = V() * sin(theta());

        // The first 6 functions are differential equations (involving yp_)
        f_[0] = -yp_[0] + omega0_ * y_[1];
        f_[1] = -yp_[1] + (1 / (2 * H_)) * ((Pm_ - D_ * y_[1]) / (1 + y_[1]) - y_[12]);
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

        return 0;
    }

    /**
     * @brief Creates the Jacobian using J and Jp.
     *
     * Values were calculated by hand and verified with MATLAB.
     *
     */
    template <class ScalarT, typename IdxT>
    int GENROU<ScalarT, IdxT>::evaluateJacobian()
    {
        J_.zeroMatrix();
        double SB = ((1 / 0.2) * (1 / 0.2)) * ((sqrt(S1_) - sqrt(S12_)) * (sqrt(S1_) - sqrt(S12_)));
        double SA = 1 - sqrt(S1_ / SB);
        const double Vr = V() * cos(theta());
        const double Vi = V() * sin(theta());

        // These variables are for ease of reading
        double df2_domega = Pm_ / ((1 + (y_[1]) * (y_[1])) * 2 * H_) - D_ / ((1 + y_[1]) * 2 * H_) + D_ * y_[1] / ((1 + (y_[1]) * (y_[1])) * 2 * H_);
        double df2_dy12 = -1 / (2 * H_);

        double df3_dId = -(Xd_ - Xdp_) * (1 - (Xdp_ - Xl_) * ((Xdp_ - Xdpp_) / ((Xdp_ - Xl_) * (Xdp_ - Xl_))));
        double df3_dksat = -y_[7];
        double df3_dpsippd = -y_[9];
        double df3_dpsipd = (Xd_ - Xdp_) * (Xdp_ - Xdpp_) / ((Xdp_ - Xl_) * (Xdp_ - Xl_));
        double df3_dEpq = -1 - (Xd_ - Xdp_) * ((Xdp_ - Xdpp_) / ((Xdp_ - Xl_) * (Xdp_ - Xl_)));

        double df6_dpsippq = y_[9] * ((Xq_ - Xl_) / (Xd_ - Xl_)) / Tq0p_;
        double df6_dEpd = -(1 / Tq0p_) * (1 + (Xq_ - Xqp_) * ((Xqp_ - Xqpp_) / ((Xqp_ - Xl_) * (Xqp_ - Xl_))));
        double df6_dpsipq = (Xq_ - Xqp_) * ((Xqp_ - Xqpp_) / ((Xqp_ - Xl_) * (Xqp_ - Xl_))) / Tq0p_;
        double df6_dksat = y_[6] * ((Xq_ - Xl_) / (Xd_ - Xl_)) / Tq0p_;
        double df6_dIq = (Xq_ - Xqp_) / Tq0p_ - (Xqp_ - Xl_) * (Xq_ - Xqp_) * ((Xqp_ - Xqpp_) / ((Xqp_ - Xl_) * (Xqp_ - Xl_))) / Tq0p_;

        double df7_dpsipq = -((Xqp_ - Xqpp_) / (Xqp_ - Xl_));
        double df7_dEpd = -((Xqpp_ - Xl_) / (Xqp_ - Xl_));

        double df8_dpsipd = (Xdp_ - Xdpp_) / (Xdp_ - Xl_);
        double df8_dEpq = (Xdpp_ - Xl_) / (Xdp_ - Xl_);

        double dTelec_dId = -Xdpp_ * y_[14] - y_[6] + y_[14] * Xdpp_;
        double dTelec_dIq = y_[7] - y_[13] * Xdpp_ + Xdpp_ * y_[13];

        // Create dF/dy
        std::vector<IdxT> rcord{0,
                                1, 1,
                                2, 2, 2, 2, 2,
                                3, 3, 3,
                                4, 4, 4,
                                5, 5, 5, 5, 5,
                                6, 6, 6,
                                7, 7, 7,
                                8, 8, 8,
                                9, 9,
                                10, 10, 10,
                                11, 11, 11,
                                12, 12, 12, 12, 12,
                                13, 13, 13,
                                14, 14, 14};

        std::vector<IdxT> ccord{1,
                                1, 12,
                                2, 3, 7, 9, 13,
                                2, 3, 13,
                                4, 5, 14,
                                4, 5, 6, 9, 14,
                                4, 5, 6,
                                2, 3, 7,
                                6, 7, 8,
                                8, 9,
                                1, 6, 10,
                                1, 7, 11,
                                6, 7, 12, 13, 14,
                                0, 13, 14,
                                0, 13, 14};

        std::vector<ScalarT> vals{omega0_,
                                  -df2_domega, df2_dy12,
                                  df3_dEpq, df3_dpsipd, df3_dpsippd, df3_dksat, df3_dId,
                                  1 / Td0pp_, -1 / Td0pp_, (Xl_ - Xdp_) / Td0pp_,
                                  -1 / Tq0pp_, 1 / Tq0pp_, (Xqp_ - Xl_) / Tq0pp_,
                                  df6_dpsipq, df6_dEpd, df6_dpsippq, df6_dksat, df6_dIq,
                                  df7_dpsipq, df7_dEpd, -1,
                                  df8_dEpq, df8_dpsipd, -1,
                                  y_[6] / sqrt((y_[7] * y_[7]) + (y_[6] * y_[6])), y_[7] / sqrt((y_[7] * y_[7]) + (y_[6] * y_[6])), -1,
                                  2.0 * SB * (y_[8] - SA), -1,
                                  -y_[6], -(1 + y_[1]), -1,
                                  y_[7], 1 + y_[1], -1,
                                  -y_[13], y_[14], -1, dTelec_dId, dTelec_dIq,
                                  Vr * cos(y_[0]) + Vi * sin(y_[0]), Rs_, -Xqpp_,
                                  -Vr * sin(y_[0]) + Vi * cos(y_[0]), Xqpp_, Rs_};
        J_.setValues(rcord, ccord, vals);

        // Create dF/dy'
        std::vector<IdxT> rcordder{0, 1, 2, 3, 4, 5};
        std::vector<IdxT> ccordder{0, 1, 2, 3, 4, 5};
        std::vector<ScalarT> valsder{-1, -1, -Td0p_, -1, -1, -1};
        COO_Matrix<ScalarT, IdxT> Jacder = COO_Matrix<ScalarT, IdxT>(rcordder, ccordder, valsder, 15, 15);

        // Perform dF/dy + \alpha dF/dy'
        J_.axpy(alpha_, Jacder);
        return 0;
    }

    /**
     * @brief This is not used, but is implemented to keep SystemModel.hpp happy.
     *
     */
    template <class ScalarT, typename IdxT>
    int GENROU<ScalarT, IdxT>::evaluateIntegrand()
    {
        return 0;
    }

    /**
     * @brief This is not used, but is implemented to keep SystemModel.hpp happy.
     *
     */
    template <class ScalarT, typename IdxT>
    int GENROU<ScalarT, IdxT>::initializeAdjoint()
    {
        return 0;
    }

    /**
     * @brief This is not used, but is implemented to keep SystemModel.hpp happy.
     *
     */
    template <class ScalarT, typename IdxT>
    int GENROU<ScalarT, IdxT>::evaluateAdjointResidual()
    {
        return 0;
    }

    /**
     * @brief This is not used, but is implemented to keep SystemModel.hpp happy.
     *
     */
    template <class ScalarT, typename IdxT>
    int GENROU<ScalarT, IdxT>::evaluateAdjointIntegrand()
    {
        return 0;
    }

    // Available template instantiations
    template class GENROU<double, long int>;
    template class GENROU<double, size_t>;

} // namespace ModelLib