// Made by Paul Moon 6/7/2024


#include <iostream>
#include <cmath>
#include <ComponentLib/PowerFlow/Bus/BaseBus.hpp>
#include <complex>
#include "GENROU.hpp"

namespace ModelLib {

/*!
 * @brief Constructor for a simple generator model
 *
 * Arguments passed to ModelEvaluatorImpl:
 * - Number of equations = 4 differential + 2 algebraic = 6
 * - Number of quadratures = 1
 * - Number of optimization parameters = 2
 */
template <class ScalarT, typename IdxT>
GENROU<ScalarT, IdxT>::GENROU(bus_type* bus, ScalarT P0, ScalarT Q0)
  : ModelEvaluatorImpl<ScalarT, IdxT>(20, 1, 2),
    H_(3.0),
    D_(0.0),
    Xq_(2.0),
    Xd_(2.1),
    Xqp_(0.5),
    Xdp_(0.3),
    Xqpp_(0.28),
    Xdpp_(0.28),
    Rs_(0.0),
    Tq0p_(0.75), 
    Td0p_(7.0), 
    Tq0pp_(0.07), 
    Td0pp_(0.073), 
    Ef_(2.9133),
    Pm_(1.0),
    X_(0.22),
    Xl_(0.13),
    omega_s_(0.0),
    omega0_(2.0*60.0*M_PI),
    S1_(0.05),
    S12_(0.2),
    omega_up_(omega_s_ + 0.0001),
    omega_lo_(omega_s_ - 0.0001),
    c_(10000.0),
    beta_(2),
    P0_(P0),
    Q0_(Q0),
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
    //std::cout << "Allocate GENROU..." << std::endl;
    tag_.resize(size_);

    return 0;
}

/**
 * @brief Initialization of the generator model
 *
 * Initialization equations are derived from example 9.2 in Power System
 * Modeling and Scripting, Federico Milano, Chapter 9, p. 225:
 * \f{eqnarray*}{
 * &~& \omega_0 = 0, \\
 * &~& \delta_0 = \tan^{-1} \left(\frac{X_q P_0 - R_s Q_0}{V_0^2 + R_s P_0 + X_q Q_0} \right) + \theta_0, \\
 * &~& \phi_0   = \delta_0 - \theta_0 + \tan^{-1} \left( \frac{Q_0}{P_0} \right), \\
 * &~& I_{d0}   = \frac{\sqrt{P_0^2 + Q_0^2}}{V_0} \sin(\phi_0), \\
 * &~& I_{q0}   = \frac{\sqrt{P_0^2 + Q_0^2}}{V_0} \cos(\phi_0), \\
 * &~& E_{d0}'  = V_0 \sin(\delta_0 - \theta_0) + R_s I_{d0} - X_q' I_{q0}, \\
 * &~& E_{q0}'  = V_0 \cos(\delta_0 - \theta_0) + R_s I_{q0} + X_d' I_{d0}
 * \f}
 *
 * The input from exciter and governor is set to the steady state value:
 * \f{eqnarray*}{
 * &~& E_{f0} = E_{q0}' + (X_d - X_d') I_{d0}, \\
 * &~& P_{m0} = E_{d0}' I_{d0} + E_{q0}' I_{q0} + ( X_q' - X_d') I_{d0} I_{q0}
 * \f}
 *
 */
template <class ScalarT, typename IdxT>
int GENROU<ScalarT, IdxT>::initialize()
{
    // std::cout << "Initialize GENROU..." << std::endl;
    //const ScalarT Vr = V()*cos(theta());
    //const ScalarT Vi = V()*sin(theta());
    //std::cout << "Initial Vr & Vi: " << Vr << " " << Vi << std::endl;

    //const ScalarT Ir = (P0_*Vr + Q0_*Vi)/((Vr*Vr)+(Vi*Vi));
    //const ScalarT Ii = -(Q0_*Vr - P0_*Vi)/((Vr*Vr)+(Vi*Vi));
    //std::cout << "Initial Ir & Ii: " << Ir << " " << Ii << std::endl;

    // Compute initial guess for the generator voltage phase
    //const ScalarT delta = atan((Xq_*P0_ - Rs_*Q0_) / (V()*V() + Rs_*P0_ + Xq_*Q0_)) + theta();

    // Compute initial guess for the generator current phase
    //const ScalarT phi   = theta() - delta - atan(Q0_/P0_);

    //std::cout << "Initial delta & phi: " << delta << " " << phi << std::endl;

    // Compute initial gueses for generator currents and potentials in d-q frame
    //const ScalarT Id = std::sqrt(P0_*P0_ + Q0_*Q0_)/V() * sin(phi);
    //const ScalarT Iq = std::sqrt(P0_*P0_ + Q0_*Q0_)/V() * cos(phi);
    //const ScalarT Id = 0.9909;
    //const ScalarT Iq = 0.3553;

    //std::cout << "Initial Id & Iq: " << Id << " " << Iq << std::endl;

   // const ScalarT Vd = V()*sin(theta() - delta) + Rs_*Id + Xqp_*Iq;
   // const ScalarT Vq = V()*cos(theta() - delta) + Rs_*Iq - Xdp_*Id;

    //std::cout << "Initial Vd & Vq: " << Vd << " " << Vq << std::endl;

    /*const ScalarT Vr = 1.0723;
    const ScalarT Vi = 0.22;
    const ScalarT Ir = 1;
    const ScalarT Ii = -0.3286;
    const ScalarT delta = 0.909;
    const ScalarT Vd = 0.7107;
    const ScalarT Vq = 0.8326;
    const ScalarT Eqp = 1.1298;
    const ScalarT Psidp = 0.9614;
    const ScalarT Psiqp = 0.6645;
    const ScalarT Edp = 0.533;
    const ScalarT Psiqpp = -0.6111;
    const ScalarT Psidpp = 1.110;
    const ScalarT Psipp = 1.2671;
    const ScalarT ksat = 1;
    const ScalarT Telec = 1;
    const ScalarT Pmech = 0.9999;*/


    //y0 = [delta_init; omega_init; E_p_q_init; psi_p_d_init; psi_p_q_init; E_p_d_init; psi_pp_q_init; psi_pp_d_init; psi_pp_init; ksat_init; Vd_init; Vq_init ; Telec_init; Id_init; Iq_init; Vr_init; Vi_init; Ir_init; Ii_init;Pmech_init]
    /*y_[0] =  delta;
    y_[1] =  omega_s_;
    y_[2] =  Eqp;
    y_[3] =  Psidp;
    y_[4] =  Psiqp;
    y_[5] =  Edp;
    y_[6] =  Psiqpp;
    y_[7] =  Psidpp;
    y_[8] =  Psipp;
    y_[9] =  ksat;
    y_[10] =  Vd;
    y_[11] =  Vq;
    y_[12] =  Telec;
    y_[13] =  Id;
    y_[14] =  Iq;
    y_[15] =  Vr;
    y_[16] =  Vi;
    y_[17] =  Ir;
    y_[18] =  Ii;
    y_[19] =  Pmech;*/
    y_[0] =  0.5273;
    y_[1] =  0;
    y_[2] =  1.1948;
    y_[3] =  1.1554;
    y_[4] =  0.2446;
    y_[5] =  0;
    y_[6] =  -0.2236;
    y_[7] =  1.1790;
    y_[8] =  1.2000;
    y_[9] =  0;
    y_[10] =  0.3494;
    y_[11] =  1.0373;
    y_[12] =  1;
    y_[13] =  0.7872;
    y_[14] =  0.6989;
    y_[15] =  1.0723;
    y_[16] =  0.22;
    y_[17] =  1.0;
    y_[18] =  -0.3286;
    y_[19] =  1.0;
    /*y_[0] =  0.9224;
    y_[1] =  0.1765;
    y_[2] =  1.2163;
    y_[3] =  1.1554;
    y_[4] =  -0.3907;
    y_[5] =  0.0;
    y_[6] =  0.3572;
    y_[7] =  0.2676;
    y_[8] =  0.4464;
    y_[9] =  2.5613;
    y_[10] =  -0.4203;
    y_[11] =  0.3149;
    y_[12] =  0.0;
    y_[13] =  0.7872;
    y_[14] =  1.0507;
    y_[15] =  0.0797;
    y_[16] =  0.2776;
    y_[17] =  1.2620;
    y_[18] =  0.3622;
    y_[19] =  1.0;*/
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
    yp_[16] = 0.0;
    yp_[17] = 0.0;
    yp_[18] = 0.0;
    yp_[19] = 0.0;

    // Initialize optimization parameters
    param_[0] = Pm_;
    param_up_[0] = 1.5;
    param_lo_[0] = 0.0;

    param_[1] = Ef_;
    param_up_[1] = 1.7;
    param_lo_[1] = 0.0;

    /*for (int i = 0; i < 20; i++) {
        std::cout << y_[i] << ", ";
    }
    std::cout << std::endl;*/

    return 0;
}

/**
 * \brief Identify differential variables.
 */
template <class ScalarT, typename IdxT>
int GENROU<ScalarT, IdxT>::tagDifferentiable()
{
    for (IdxT i=6; i < size_; ++i)
    {
        tag_[i] = false;
    }
    tag_[0] = true;
    tag_[1] = true;
    tag_[2] = true;
    tag_[3] = true;
    tag_[4] = true;
    tag_[5] = true;

    

    return 0;
}

/**
 * @brief Computes residual vector for the generator model.
 *
 * Residual equations are given per model in Power System Modeling and
 * Scripting, Federico Milano, Chapter 15, p. 334:
 * \f{eqnarray*}{
 * f_0: &~& \dot{\delta} -\omega_b (\omega - \omega_s), \\
 * f_1: &~& 2H/\omega_s \dot{\omega} - L_m(P_m) + E_q' I_q + E_d' I_d + (X_q' - X_d')I_d I_q  + D (\omega - \omega_s), \\
 * f_2: &~& T_{q0}' \dot{E}_d' + E_d' - (X_q - X_q')I_q, \\
 * f_3: &~& T_{d0}' \dot{E}_q' + E_q' + (X_d - X_d')I_d - E_f, \\
 * f_4: &~& R_s I_d - X_q' I_q + V \sin(\delta - \theta) - E_d', \\
 * f_5: &~& R_s I_q + X_d' I_d + V \cos(\delta - \theta) - E_q',
 * \f}
 * where \f$ \Omega_b \f$ is the synchronous frequency in [rad/s], and
 * overdot denotes time derivative.
 *
 * Generator injection active and reactive power are
 * \f{eqnarray*}{
 * P_g &=& E_d' I_d + E_q' I_q + (X_q' - X_d') I_d I_q - R_s (I_d^2 + I_q^2), \\
 * Q_q &=& E_q' I_d - E_d' I_q - X_q' I_q^2 - X_d' I_d^2, \\
 * \f}
 * respectively.
 *
 * State variables are:
 * \f$ y_0 = \omega \f$, \f$ y_1 = \delta \f$, \f$ y_2 = E_d' \f$, \f$ y_3 = E_q' \f$,
 * \f$ y_4 = I_d \f$, \f$ y_5 = I_q \f$.
 *
 */
template <class ScalarT, typename IdxT>
int GENROU<ScalarT, IdxT>::evaluateResidual()
{
    double SA=std::log(S12_/S1_)/std::log(1.2);
    double SB=S1_;
    std::complex<double> numerator(y_[15], y_[16]);   // numerator is y(16) + 1i * y(17)
    std::complex<double> denominator(0.0, X_);   // denominator is 1i * X (1i is imaginary unit in MATLAB)
    std::complex<double> division_result = numerator / denominator;

    double real_part = division_result.real();  // real part of (numerator / denominator)
    double imag_part = division_result.imag();  // imaginary part of (numerator / denominator)


    //std::cout << "SA: " << SA << " SB: " << SB << std::endl;
    // std::cout << "Evaluate residual for GENROU..." << std::endl;
    f_[0] = -yp_[0] + omega0_*y_[1]; //f1 = -yp(1)+y(2)*w0_1;
    f_[1] = -yp_[1] + (1/(2*H_))*((Pm_ - D_*y_[1])/(1 + y_[1]) - y_[12]);  //f2 = -yp(2) + (1/(2*H_1))*((Pmech_init - D_1*y(2))/(1 + y(2)) - y(13));
    f_[2] = -yp_[2]*Td0p_ + Ef_ - (y_[2] + (Xd_-Xdp_)*y_[13]+((Xdp_-Xdpp_)/((Xdp_-Xl_)*(Xdp_-Xl_)))*(y_[2]-y_[3]-(Xdp_-Xl_)*y_[13])) + y_[7]*y_[9]; //f3 = -yp(3)*T_p_d0_1 + Efd - (y(3) + (Xd_1-X_p_d_1)*(y(14) + ((X_p_d_1-X_pp_d_1)/(X_p_d_1 - Xl_1)^2)*(y(3)-y(4)-(X_p_d_1-Xl_1)*y(14))) + y(8)*y(10));
    f_[3] = -yp_[3]+(1/Td0pp_)*(y_[2]-y_[3]-(Xdp_-Xl_)*y_[13]); //f4 = -yp(4)+(1/T_pp_d0_1)*(y(3)-y(4)-(X_p_d_1-Xl_1)*y(14));
    f_[4] = -yp_[4]+(1/Tq0pp_)*(y_[5]-y_[4]+(Xqp_-Xl_)*y_[14]); //f5 = -yp(5)+(1/T_pp_q0_1)*(y(6)-y(5)+(X_p_q_1-Xl_1)*y(15));
    f_[5] = -yp_[5] + (1/Tq0p_)*(-y_[5]+((Xq_-Xl_)/(Xd_-Xl_)*y_[6]*y_[9] + (Xq_-Xqp_)*(y_[14]-((Xqp_-Xqpp_)/((Xqp_-Xl_)*(Xqp_-Xl_))*y_[5]+y_[14])*(Xqp_-Xl_)-y_[4]))); //f6 = -yp(6) + (1/T_p_q0_1)*(-y(6)+((Xq_1-Xl_1)/(Xd_1-Xl_1))*y(7)*y(10) +(Xq_1-X_p_q_1)*(y(15)-((X_p_q_1-X_pp_q_1)/(X_p_q_1-Xl_1)^2)*(y(6)+y(15)*(X_p_q_1-Xl_1)-y(5))));

    f_[6] = -y_[6] - y_[4]*((Xqp_-Xqpp_)/(Xqp_-Xl_))-y_[5]*((Xqpp_-Xl_)/(Xqp_-Xl_)); //f7 = -y(7) - y(5)*((X_p_q_1-X_pp_q_1)/(X_p_q_1 - Xl_1))- y(6)*((X_pp_q_1 - Xl_1)/(X_p_q_1 - Xl_1));
    f_[7] = -y_[7] - y_[3]*((Xdp_-Xdpp_)/(Xdp_-Xl_))+y_[2]*((Xdpp_-Xl_)/(Xdp_-Xl_)); //f8 = -y(8) - y(4)*((X_p_d_1-X_pp_d_1)/(X_p_d_1 - Xl_1)) + y(3)*((X_pp_d_1 - Xl_1)/(X_p_d_1 - Xl_1));
    f_[8] = -y_[8] + sqrt((y_[7]*y_[7])+(y_[6]*y_[6]));//f9 = -y(9) + sqrt(y(8)^2+y(7)^2);
    f_[9] = -y_[9] + SB*((y_[8]-SA)*(y_[8]-SA)); //f10 = -y(10) + SB * (y(9) - SA)^2;
    f_[10] = -y_[10] - y_[6] * (1+y_[1]); //f11 = -y(11) - y(7) * (1+y(2));
    f_[11] = -y_[11] + y_[7] * (1+y_[1]);  //f12 = -y(12) + y(8) * (1+y(2));
    f_[12] = -y_[12] + y_[14] * (y_[7] - y_[13] * Xdpp_) - y_[13] * (y_[6] - y_[14] * Xdpp_); //f13 = -y(13) + y(15) * (y(8) - y(14) * X_pp_d_1) - y(14) * (y(7) - y(15) * X_pp_d_1);
    f_[13] = -y_[13] + y_[17] * sin(y_[0]) - y_[18] * cos(y_[0]); //f14 = -y(14) + y(18) * sin(y(1)) - y(19) * cos(y(1)); % prev version: -y(14) + Ireal * sin(y(1)) - Iimag * cos(y(1)); %update Ireal and Iimag along with state variables each iteration
    f_[14] = -y_[14] + y_[17] * cos(y_[0]) + y_[18] * sin(y_[0]); //f15 = -y(15) + y(18) * cos(y(1)) + y(19) * sin(y(1)); % prev: -y(15) + Ireal * cos(y(1)) + Iimag * sin(y(1));
    f_[15] = -y_[10] + y_[15]*sin(y_[0]) - y_[16]*cos(y_[0])+y_[13]*Rs_ - y_[14]*Xqpp_; //f16 = -y(11) + y(16)*sin(y(1)) - y(17)*cos(y(1))+y(14)*Ra_1 - y(15)*X_pp_q_1;
    f_[16] = -y_[11] + y_[15]*cos(y_[0]) + y_[16]*sin(y_[0]) + y_[13]*Xqpp_ + y_[14]*Rs_;  //f17 = -y(12) + y(16)*cos(y(1)) + y(17)*sin(y(1)) + y(14)*X_pp_q_1 + y(15)*Ra_1;
    f_[17] = -y_[17] + real_part;// y_[16]/X_; //f18 = -y(18)+real((y(16) + 1i * y(17))/(1i*X));
    f_[18] = -y_[18] + imag_part;//y_[15]/X_; //f19 = -y(19)+imag((y(16) + 1i * y(17))/(1i*X));
    f_[19] = -Pm_ + y_[19]; //f20 = -Pmech_init + y(20);

    return 0;
}

template <class ScalarT, typename IdxT>
int GENROU<ScalarT, IdxT>::evaluateJacobian()
{
    //std::cout << "Evaluate Jacobian for GENROU..." << std::endl;
    J_.zeroMatrix();
    double SA=log(S12_/S1_)/log(1.2);
    double SB=S1_;

    // These variables are for ease of reading
    double df2_domega= y_[19]/((1+(y_[1])*(y_[1]))*2*H_) - D_/((1+y_[1])*2*H_) + D_*y_[1]/((1+(y_[1])*(y_[1]))*2*H_);

    double df3_dId = -(Xd_ - Xdp_) * (1 - (Xdp_ - Xl_) * ((Xdp_ - Xdpp_) / ((Xdp_ - Xl_)*(Xdp_ - Xl_))));
    double df3_dksat = -y_[7];
    double df3_dpsippd = -y_[9];
    double df3_dpsipd = (Xd_ - Xdp_)*(Xdp_ - Xdpp_)/((Xdp_ - Xl_)*(Xdp_ - Xl_));
    double df3_dEpq = -1 - (Xd_ - Xdp_) * ((Xdp_ - Xdpp_)/((Xdp_ - Xl_)*(Xdp_ - Xl_)));

    double df6_dpsippq = y_[9]*((Xq_-Xl_)/(Xd_-Xl_))/Tq0p_;
    double df6_dEpd = -(1/Tq0p_) * (1 + (Xq_-Xqp_)*((Xqp_ - Xqpp_)/((Xqp_ - Xl_)*(Xqp_ - Xl_))));
    double df6_dpsipq = (Xq_-Xqp_)*((Xqp_ - Xqpp_)/((Xqp_ - Xl_)*(Xqp_ - Xl_)))/Tq0p_;
    double df6_dksat = y_[6]*((Xq_-Xl_)/(Xd_-Xl_))/Tq0p_;
    double df6_dIq = (Xq_-Xqp_)/Tq0p_ - (Xqp_-Xl_)*(Xq_-Xqp_)*((Xqp_ - Xqpp_)/((Xqp_ - Xl_)*(Xqp_ - Xl_)))/Tq0p_;

    double df7_dpsipq = -((Xqp_ - Xqpp_) / (Xqp_ - Xl_));
    double df7_dEpd = - ((Xqpp_ - Xl_)/(Xqp_ - Xl_));

    double df8_dpsipd = -(Xdp_-Xdpp_)/(Xdp_ - Xl_);
    double df8_dEpq = (Xdpp_ - Xl_)/(Xdp_ - Xl_);

    double dTelec_dId = -Xdpp_*y_[14] - y_[6] + y_[14]*Xdpp_;
    double dTelec_dIq = y_[7] - y_[13]*Xdpp_ + Xdpp_*y_[13];

    //Create dF/dy
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
    13, 13, 13, 13,
    14, 14, 14, 14,
    15, 15, 15, 15, 15, 15,
    16, 16, 16, 16, 16, 16,
    17, 17,
    18, 18,
    19};
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
    0, 13, 17, 18,
    0, 14, 17, 18,
    0, 10, 13, 14, 15, 16,
    0, 11, 13, 14, 15, 16,
    16, 17,
    15, 18,
    19};
    std::vector<ScalarT> vals{omega0_,
    -df2_domega, (-1/(2*H_)),
    df3_dEpq, df3_dpsipd, df3_dpsippd, df3_dksat, df3_dId,
    1/Td0pp_, -1/Td0pp_, (Xl_-Xdp_)/Td0pp_,
    -1/Tq0pp_, 1/Tq0pp_, (Xqp_-Xl_)/Tq0pp_,
    df6_dpsipq, df6_dEpd, df6_dpsippq, df6_dksat, df6_dIq,
    df7_dpsipq, df7_dEpd, -1,
    df8_dEpq, df8_dpsipd, -1,
    y_[6]/sqrt((y_[7]*y_[7])+(y_[6]*y_[6])), y_[7]/sqrt((y_[7]*y_[7])+(y_[6]*y_[6])), -1,
    2*SB*(y_[8] - SA), -1,
    -y_[6], -(1+y_[1]), -1,
    y_[7], 1+y_[1], -1,
    -y_[13], y_[14], -1, dTelec_dId, dTelec_dIq,
    y_[17]*cos(y_[0])+y_[18]*sin(y_[0]), -1, sin(y_[0]), -cos(y_[0]),
    -y_[17]*sin(y_[0])+y_[18]*cos(y_[0]), -1, cos(y_[0]), sin(y_[0]),
    y_[15]*cos(y_[0])+y_[16]*sin(y_[0]), -1, Rs_, -Xqpp_, sin(y_[0]), -cos(y_[0]),
    -y_[15]*sin(y_[0])+y_[16]*cos(y_[0]), -1, Xqpp_, Rs_, cos(y_[0]), sin(y_[0]),
    1.0/X_, -1,
    -1.0/X_, -1,
    1};
    J_.setValues(rcord, ccord, vals);

    //Create dF/dy'
    std::vector<IdxT> rcordder{0,1,2,3,4,5};
    std::vector<IdxT> ccordder{0,1,2,3,4,5};
    std::vector<ScalarT> valsder{-1,-1,-Td0p_,-1,-1,-1};
    COO_Matrix<ScalarT,IdxT> Jacder = COO_Matrix<ScalarT, IdxT>(rcordder, ccordder, valsder,20,20);
    
    //Perform dF/dy + \alpha dF/dy'
    J_.axpy(alpha_, Jacder);
    return 0;
}

template <class ScalarT, typename IdxT>
int GENROU<ScalarT, IdxT>::evaluateIntegrand()
{
    // std::cout << "Evaluate Integrand for GENROU..." << std::endl;
    g_[0] = frequencyPenalty(y_[1]);
    return 0;
}

template <class ScalarT, typename IdxT>
int GENROU<ScalarT, IdxT>::initializeAdjoint()
{
    //std::cout << "Initialize adjoint for GENROU..." << std::endl;
    for (IdxT i=0; i<size_; ++i)
    {
        yB_[i] = 0.0;
        ypB_[i] = 0.0;
    }
    ypB_[1] = frequencyPenaltyDer(y_[1]);

    return 0;
}


/**
 * @brief Computes adjoint residual vector for the generator model.
 *
 * Adjoint residual equations are given as:
 * \f{eqnarray*}{
 * f_{B0}: &~& \dot{y}_{B0} - y_{B4} V \cos(\delta - \theta) + y_{B5} V \sin(\delta - \theta), \\
 * f_{B1}: &~& 2H/\omega_s \dot{y}_{B1} + y_{B0} \omega_b - y_{B1} D + y_{B9} (1 - T_2/T_1) - y_{B10} K T_2/T_1 + g_{\omega}(\omega), \\
 * f_{B2}: &~& T_{q0}' \dot{y}_{B2} - y_{B1} I_d - y_{B2} + y_{B4} + y_{B6} I_d - y_{B7} I_q, \\
 * f_{B3}: &~& T_{d0}' \dot{y}_{B3} - y_{B1} I_q - y_{B3} + y_{B5} + y_{B6} I_q + y_{B7} I_d, \\
 * f_{B4}: &~& -y_{B1} (E_d' + (-X_d'+X_q') I_q) - y_{B3} (X_d - X_d') - y_{B4} R_s - y_{B5} X_d' + y_{B6} (E_d' + (X_q' - X_d') I_q - 2 R_s I_d) + y_{B7} (E_q' - 2 X_d' I_d), \\
 * f_{B5}: &~& -y_{B1} (E_q' + (-X_d'+X_q') I_d) + y_{B2} (X_q - X_q') + y_{B4} X_q' - y_{B5} R_s + y_{B6} (E_q' + (X_q' - X_d') I_d - 2 R_s I_q) - y_{B7} (E_d' + 2 X_q' I_q). \\
 * \f}
 *
 */
template <class ScalarT, typename IdxT>
int GENROU<ScalarT, IdxT>::evaluateAdjointResidual()
{
    // std::cout << "Evaluate adjoint residual for GENROU..." << std::endl;
    ScalarT sinPhi = sin(delta() - theta());
    ScalarT cosPhi = cos(delta() - theta());

    // Generator adjoint
    fB_[0] = ypB_[0] - yB_[4]*V()*cosPhi + yB_[5]*V()*sinPhi;
    fB_[1] = 2.0*H_/omega_s_*ypB_[1] + yB_[0]*omega0_ - yB_[1]*D_ + frequencyPenaltyDer(omega());
    fB_[2] = Tq0p_*ypB_[2] - yB_[1]*Id() - yB_[2] + yB_[4];
    fB_[3] = Td0p_*ypB_[3] - yB_[1]*Iq() - yB_[3] + yB_[5];
    fB_[4] = -yB_[1]*(Edp() + (Xqp_ - Xdp_)*Iq()) - yB_[3]*(Xd_ - Xdp_) - yB_[4]*Rs_ - yB_[5]*Xdp_;
    fB_[5] = -yB_[1]*(Eqp() + (Xqp_ - Xdp_)*Id()) + yB_[2]*(Xq_ - Xqp_) + yB_[4]*Xqp_ - yB_[5]*Rs_;

    return 0;
}

// template <class ScalarT, typename IdxT>
// int GENROU<ScalarT, IdxT>::evaluateAdjointJacobian()
// {
//     std::cout << "Evaluate adjoint Jacobian for GENROU..." << std::endl;
//     std::cout << "Adjoint Jacobian evaluation not implemented!" << std::endl;
//     return 0;
// }

template <class ScalarT, typename IdxT>
int GENROU<ScalarT, IdxT>::evaluateAdjointIntegrand()
{
    // std::cout << "Evaluate adjoint Integrand for GENROU..." << std::endl;
    gB_[0] = yB_[1];
    gB_[1] = yB_[3];

    return 0;
}


//
// Private functions
//

/**
 * Generator active power Pg.
 *
 * \f[ P_g = E_q' I_q + E_d' I_d + (X_q' - X_d') I_q I_d - R_a (I_d^2 + I_q^2) \f]
 *
 */
template <class ScalarT, typename IdxT>
ScalarT GENROU<ScalarT, IdxT>::Pg()
{
    return 1;
}

/**
 * Generator reactive power Qg.
 *
 * \f[ Q_g = E_q' I_d - E_d' I_q - X_d' I_d^2 - X_q' I_q^2 \f]
 */
template <class ScalarT, typename IdxT>
ScalarT GENROU<ScalarT, IdxT>::Qg()
{
    return 1;
}

/**
 * Frequency penalty is used as the objective function for the generator model.
 */
template <class ScalarT, typename IdxT>
ScalarT GENROU<ScalarT, IdxT>::frequencyPenalty(ScalarT omega)
{
    return 0;
}

/**
 * Derivative of frequency penalty cannot be written in terms of min/max functions.
 * Need to expand conditional statements instead.
 */
template <class ScalarT, typename IdxT>
ScalarT GENROU<ScalarT, IdxT>::frequencyPenaltyDer(ScalarT omega)
{
    return 0.0;
}



// Available template instantiations
template class GENROU<double, long int>;
template class GENROU<double, size_t>;


} // namespace ModelLib