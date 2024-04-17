


#include <iostream>
#include <cmath>
#include <vector>
#include "DistributedGenerator.hpp"

namespace ModelLib {


/*!
 * @brief Constructor for a Discrete Generator
 * @todo Maybe have parameters be templated in. Variables cannot be changed and are unlikely to. Allows for compile time optimizations
 *
 * Calls default ModelEvaluatorImpl constructor.
 */

template <class ScalarT, typename IdxT>
DistributedGenerator<ScalarT, IdxT>::DistributedGenerator(IdxT id, DistributedGeneratorParameters<ScalarT,IdxT> parm, bool reference_frame)
  : wb_(parm.wb_),
    wc_(parm.wc_),
    mp_(parm.mp_),
    Vn_(parm.Vn_), 
    nq_(parm.nq_), 
    F_(parm.F_), 
    Kiv_(parm.Kiv_), 
    Kpv_(parm.Kpv_), 
    Kic_(parm.Kic_), 
    Kpc_(parm.Kpc_), 
    Cf_(parm.Cf_), 
    rLf_(parm.rLf_), 
    Lf_(parm.Lf_), 
    rLc_(parm.rLc_), 
    Lc_(parm.Lc_), 
    refframe_(reference_frame)
{
    // internals [\delta_i, Pi, Qi, phi_di, phi_qi, gamma_di, gamma_qi, il_di, il_qi, vo_di, vo_qi, io_di, io_qi]
    // externals [\omega_ref, vba_out, vbb_out]
    size_ = 16;
    n_intern_ = 13;
    n_extern_ = 3;
    extern_indices_ = {0,1,2};
    idc_ = id;
}

template <class ScalarT, typename IdxT>
DistributedGenerator<ScalarT, IdxT>::~DistributedGenerator()
{
}

/*!
 * @brief allocate method computes sparsity pattern of the Jacobian.
 */
template <class ScalarT, typename IdxT>
int DistributedGenerator<ScalarT, IdxT>::allocate()
{
    y_.resize(size_);
    yp_.resize(size_);
    f_.resize(size_);
    
    return 0;
}

/**
 * Initialization of the grid model
 */
template <class ScalarT, typename IdxT>
int DistributedGenerator<ScalarT, IdxT>::initialize()
{
    return 0;
}

/*
 * \brief Identify differential variables
 */
template <class ScalarT, typename IdxT>
int DistributedGenerator<ScalarT, IdxT>::tagDifferentiable()
{
    return 0;
}

/**
 * @brief Contributes to the resisdual of the Distributed Generator
 *
 */
template <class ScalarT, typename IdxT>
int DistributedGenerator<ScalarT, IdxT>::evaluateResidual()
{
    //   ### Externals Componenets ###

    ScalarT omega = wb_ - mp_ * y_[4];
    //ref common ref motor angle
    /// @todo fix boolian conditional, unclear result
    if (refframe_)
    {
        f_[0] = omega - y_[0];
    }
    else
    {
        f_[0] = 0.0;
    }
    

    //output
    //current transformed to common frame
    f_[1] = cos(y_[3]) * y_[14] - sin(y_[3]) * y_[15];
    f_[2] = sin(y_[3]) * y_[14] + cos(y_[3]) * y_[15];

    //Take incoming voltages to current rotator reference frame
    ScalarT vbd_in = cos(y_[3]) * y_[1] + sin(y_[3]) * y_[2];
    ScalarT vbq_in = -sin(y_[3]) * y_[1] + cos(y_[3]) * y_[2];

    // ### Internal Componenets ## 
    // Rotor difference angle
    f_[3] = -yp_[3] + omega - y_[0];

    // P and Q equations
    f_[4] = -yp_[4] + wc_ * (y_[12] * y_[14] + y_[13] * y_[15] - y_[4]);
    f_[5] = -yp_[5] + wc_ * (-y_[12] * y_[15] + y_[13] * y_[14] - y_[5]);

    //Voltage control
    ScalarT vod_star = Vn_ - nq_ * y_[5];
    ScalarT voq_star = 0.0;

    f_[6] = -yp_[6] + vod_star - y_[12];
    f_[7] = -yp_[7] + voq_star - y_[13];


    ScalarT ild_star = F_ * y_[14] - wb_ * Cf_ * y_[13] + Kpv_ * (vod_star - y_[12]) + Kiv_ * y_[6];
    ScalarT ilq_star = F_ * y_[15] + wb_ * Cf_ * y_[12] + Kpv_ * (voq_star - y_[13]) + Kiv_ * y_[7];

    //Current control
    f_[8] = -yp_[8] + ild_star - y_[10];
    f_[9] = -yp_[9] + ilq_star - y_[11];

    ScalarT vid_star = -wb_ * Lf_ * y_[11] + Kpc_ * (ild_star - y_[10]) + Kic_ * y_[8];
    ScalarT viq_star = wb_ * Lf_ * y_[10] + Kpc_ * (ilq_star - y_[11]) + Kic_ * y_[9]; 
    
    //Output LC Filter 
    f_[10] = -yp_[10] - (rLf_ / Lf_) * y_[10] + omega * y_[11] +  (vid_star - y_[12]) / Lf_;
    f_[11] = -yp_[11] - (rLf_ / Lf_) * y_[11] - omega * y_[10] +  (viq_star - y_[13]) / Lf_;

    f_[12] = -yp_[12] + omega * y_[13] + (y_[10] - y_[14]) / Cf_;
    f_[13] = -yp_[13] - omega * y_[12] + (y_[11] - y_[15]) / Cf_;

    //Output Connector
    f_[14] = -yp_[14] - (rLc_ / Lc_) * y_[14] + omega * y_[15] + (y_[12] - vbd_in) / Lc_;
    f_[15] = -yp_[15] - (rLc_ / Lc_) * y_[15] - omega * y_[14] + (y_[13] - vbq_in) / Lc_;
    return 0;
}

/**
 * @brief Compute the jacobian of the DistributedGenerator for iteration. dF/dy - \alpha dF/dy'
 * 
 * The matrix dF/dy should be
 *  
[ 0,           0,           0,                             0,       0,                0,            0,            0,      0,      0,                 0,                 0,                 0,                 0,          0,          0]
[ 0,           0,           0,                             0,       0,                0,            0,            0,      0,      0,                 0,                 0,                 0,                 0,          0,          0]
[ 0,           0,           0,                             0,       0,                0,            0,            0,      0,      0,                 0,                 0,                 0,                 0,          0,          0]
[-1,           0,           0,                             0,     -mp,                0,            0,            0,      0,      0,                 0,                 0,                 0,                 0,          0,          0]
[ 0,           0,           0,                             0,     -wc,                0,            0,            0,      0,      0,                 0,                 0,            wc*x15,            wc*x16,     wc*x13,     wc*x14]
[ 0,           0,           0,                             0,       0,              -wc,            0,            0,      0,      0,                 0,                 0,           -wc*x16,            wc*x15,     wc*x14,    -wc*x13]
[ 0,           0,           0,                             0,       0,              -nq,            0,            0,      0,      0,                 0,                 0,                -1,                 0,          0,          0]
[ 0,           0,           0,                             0,       0,                0,            0,            0,      0,      0,                 0,                 0,                 0,                -1,          0,          0]
[ 0,           0,           0,                             0,       0,          -Kpv*nq,          Kiv,            0,      0,      0,                -1,                 0,              -Kpv,            -Cf*wb,          F,          0]
[ 0,           0,           0,                             0,       0,                0,            0,          Kiv,      0,      0,                 0,                -1,             Cf*wb,              -Kpv,          0,          F]
[ 0,           0,           0,                             0, -mp*x12, -(Kpc*Kpv*nq)/Lf, (Kiv*Kpc)/Lf,            0, Kic/Lf,      0, - Kpc/Lf - rLf/Lf,            -mp*x5, -(Kpc*Kpv + 1)/Lf,   -(Cf*Kpc*wb)/Lf, (F*Kpc)/Lf,          0]
[ 0,           0,           0,                             0,  mp*x11,                0,            0, (Kiv*Kpc)/Lf,      0, Kic/Lf,             mp*x5, - Kpc/Lf - rLf/Lf,    (Cf*Kpc*wb)/Lf, -(Kpc*Kpv + 1)/Lf,          0, (F*Kpc)/Lf]
[ 0,           0,           0,                             0, -mp*x14,                0,            0,            0,      0,      0,              1/Cf,                 0,                 0,        wb - mp*x5,      -1/Cf,          0]
[ 0,           0,           0,                             0,  mp*x13,                0,            0,            0,      0,      0,                 0,              1/Cf,        mp*x5 - wb,                 0,          0,      -1/Cf]
[ 0, -cos(x4)/Lc, -sin(x4)/Lc, -(x3*cos(x4) - x2*sin(x4))/Lc, -mp*x16,                0,            0,            0,      0,      0,                 0,                 0,              1/Lc,                 0,    -rLc/Lc, wb - mp*x5]
[ 0,  sin(x4)/Lc, -cos(x4)/Lc,  (x2*cos(x4) + x3*sin(x4))/Lc,  mp*x15,                0,            0,            0,      0,      0,                 0,                 0,                 0,              1/Lc, mp*x5 - wb,    -rLc/Lc]
 * 'Generated from MATLAB symbolic'
 * 
 * @tparam ScalarT 
 * @tparam IdxT 
 * @return int 
 */
template <class ScalarT, typename IdxT>
int DistributedGenerator<ScalarT, IdxT>::evaluateJacobian()
{
    J_.zeroMatrix();
    //Create dF/dy'
    std::vector<IdxT> rcordder(13);
    std::vector<ScalarT> valsder(13, -1.0);
    for (int i = 0; i < 13; i++)
    {
        rcordder[i] = i + 3;
    }
    COO_Matrix<ScalarT,IdxT> Jacder = COO_Matrix<ScalarT, IdxT>(rcordder, rcordder, valsder,16,16);

    std::vector<IdxT> ctemp{};
    std::vector<IdxT> rtemp{};
    std::vector<ScalarT> valtemp{};



    //Create dF/dy
    //r = 1

    ctemp = {3, 14, 15};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++) rtemp.push_back(1);
    valtemp = { - sin(y_[3]) * y_[14] - cos(y_[3]) * y_[15], cos(y_[3]),-sin(y_[3])};
    J_.setValues(rtemp, ctemp, valtemp);

    //r = 2

    ctemp = {3, 14, 15};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++) rtemp.push_back(2);
    valtemp = { cos(y_[3]) * y_[14] - sin(y_[3]) * y_[15], sin(y_[3]),cos(y_[3])};
    J_.setValues(rtemp, ctemp, valtemp);

    //r = 3
    
    ctemp = {0, 4};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++) rtemp.push_back(3);
    valtemp = {-1.0, -mp_};
    J_.setValues(rtemp, ctemp, valtemp);
    
    //r = 0
    if (refframe_)
    {
        ctemp = {0, 4};
        rtemp.clear();
        for (size_t i = 0; i < ctemp.size(); i++) rtemp.push_back(0);
        valtemp = {-1.0, -mp_};
        J_.setValues(rtemp, ctemp, valtemp);
    }
    

    //r = 4
    ctemp = {4, 12, 13, 14, 15};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++) rtemp.push_back(4);
    valtemp = {-wc_, wc_*y_[14], wc_*y_[15], wc_*y_[12], wc_*y_[13]};
    J_.setValues(rtemp, ctemp, valtemp);

    //r = 5
    ctemp = {5, 12, 13, 14, 15};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++) rtemp.push_back(5);
    valtemp = {-wc_, -wc_*y_[15], wc_*y_[14], wc_*y_[13], -wc_*y_[12]};
    J_.setValues(rtemp, ctemp, valtemp);

    //r = 6
    ctemp = {5, 12};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++) rtemp.push_back(6);
    valtemp = {-nq_, -1.0};
    J_.setValues(rtemp, ctemp, valtemp);
    
    //r = 7
    ctemp = {13};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++) rtemp.push_back(7);
    valtemp = {-1.0};
    J_.setValues(rtemp, ctemp, valtemp);

    //r = 8
    ctemp = {5,6,10,12,13,14};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++) rtemp.push_back(8);
    valtemp = {-Kpv_*nq_, Kiv_, -1.0, -Kpv_, -Cf_*wb_, F_};
    J_.setValues(rtemp, ctemp, valtemp);

    //r = 9
    ctemp = {7, 11, 12, 13, 15};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++) rtemp.push_back(9);
    valtemp = {Kiv_, -1.0, Cf_*wb_,-Kpv_,F_};
    J_.setValues(rtemp, ctemp, valtemp);
    
    //r = 10
    ctemp = {4, 5, 6, 8, 10, 11, 12, 13, 14};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++) rtemp.push_back(10);
    valtemp = {-mp_ * y_[11], -(Kpc_ * Kpv_ * nq_) / Lf_, (Kpc_ * Kiv_) / Lf_, Kic_ / Lf_, -(Kpc_ + rLf_) / Lf_, -mp_ * y_[4], -(Kpc_ * Kpv_ + 1.0) / Lf_, -(Cf_ * Kpc_ * wb_) / Lf_, (F_ * Kpc_) / Lf_};
    J_.setValues(rtemp, ctemp, valtemp);

    //r = 11
    ctemp = {4, 7, 9, 10, 11, 12, 13, 15};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++) rtemp.push_back(11);
    valtemp = {mp_ * y_[10], (Kiv_ * Kpc_) / Lf_, Kic_ / Lf_, mp_ * y_[4], -(Kpc_ + rLf_) / Lf_, (Cf_ * Kpc_ * wb_) / Lf_, -(Kpc_ * Kpv_ + 1.0) / Lf_, (F_ * Kpc_) / Lf_};
    J_.setValues(rtemp, ctemp, valtemp);

    //r = 12
    ctemp = {4, 10, 13, 14};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++) rtemp.push_back(12);
    valtemp = {-mp_ * y_[13], 1.0 / Cf_, wb_ - mp_ * y_[4], -1.0 / Cf_};
    J_.setValues(rtemp, ctemp, valtemp);

    
    //r = 13
    ctemp = {4, 11, 12, 15};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++) rtemp.push_back(13);
    valtemp = {mp_ * y_[12], 1.0 / Cf_, -wb_ + mp_ * y_[4], -1.0 / Cf_};
    J_.setValues(rtemp, ctemp, valtemp);

    
    //r = 14
    ctemp = {1, 2, 3, 4, 12, 14, 15};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++) rtemp.push_back(14);
    valtemp = {(1.0/Lc_) * -cos(y_[3]) , (1.0/Lc_) * -sin(y_[3]) , (1.0/Lc_) * (sin(y_[3]) * y_[1] - cos(y_[3]) * y_[2]), -mp_ * y_[15], 1.0 / Lc_, -rLc_ / Lc_, wb_ - mp_ * y_[4]};
    J_.setValues(rtemp, ctemp, valtemp);

    
    //r = 15
    ctemp = {1, 2, 3, 4, 13, 14, 15};
    rtemp.clear();
    for (size_t i = 0; i < ctemp.size(); i++) rtemp.push_back(15);
    valtemp = {(1.0/Lc_) * sin(y_[3]) , (1.0/Lc_) * -cos(y_[3]), (1.0/Lc_) * (cos(y_[3]) * y_[1] + sin(y_[3]) * y_[2]), mp_ * y_[14], 1.0 / Lc_, -wb_ + mp_ * y_[4], -rLc_ / Lc_};
    J_.setValues(rtemp, ctemp, valtemp);


    //Perform dF/dy + \alpha dF/dy'

    J_.axpy(alpha_, Jacder);

    return 0;
}

template <class ScalarT, typename IdxT>
int DistributedGenerator<ScalarT, IdxT>::evaluateIntegrand()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int DistributedGenerator<ScalarT, IdxT>::initializeAdjoint()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int DistributedGenerator<ScalarT, IdxT>::evaluateAdjointResidual()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int DistributedGenerator<ScalarT, IdxT>::evaluateAdjointIntegrand()
{
    return 0;
}




// Available template instantiations
template class DistributedGenerator<double, long int>;
template class DistributedGenerator<double, size_t>;


} //namespace ModelLib

