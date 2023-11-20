


#include <iostream>
#include <cmath>
#include <vector>
#include "DiscreteGenerator.hpp"

namespace ModelLib {


/*!
 * @brief Constructor for a Discrete Generator
 * @todo Maybe have parameters be templated in. Variables cannot be changed and are unlikely to. Allows for compile time optimizations
 *
 * Calls default ModelEvaluatorImpl constructor.
 */

template <class ScalarT, typename IdxT>
DiscreteGenerator<ScalarT, IdxT>::DiscreteGenerator(IdxT id, DiscreteGeneratorParameters<ScalarT,IdxT> parm)
  :  wb_(parm.wb), wc_(parm.wc), mp_(parm.mp), refmp_(parm.refmp), Vn_(parm.Vn), nq_(parm.nq), F_(parm.F), Kiv_(parm.Kiv), Kpv_(parm.Kpv), Kic_(parm.Kic), Kpc_(parm.Kpc), Cf_(parm.Cf), rLf_(parm.rLf), Lf_(parm.Lf), rLc_(parm.rLc), Lc_(parm.Lc) 
{
    // internals [delta_i, Pi, Qi, phi_di, phi_qi, gamma_di, gamma_qi, il_di, il_qi, vo_di, vo_qi, io_di, io_qi]
    // externals [Pref, Pbus, QBus]
	this->size_ = 16;
    this->n_intern = 13;
    this->n_extern = 3;
    this->extern_indices = {0,1,2};
    this->idc_ = id;
}

template <class ScalarT, typename IdxT>
DiscreteGenerator<ScalarT, IdxT>::~DiscreteGenerator()
{
}

/*!
 * @brief allocate method computes sparsity pattern of the Jacobian.
 */
template <class ScalarT, typename IdxT>
int DiscreteGenerator<ScalarT, IdxT>::allocate()
{
    this->y_.resize(this->size_);
    this->yp_.resize(this->size_);
    this->f_.resize(this->size_);
    
    return 0;
}

/**
 * Initialization of the grid model
 */
template <class ScalarT, typename IdxT>
int DiscreteGenerator<ScalarT, IdxT>::initialize()
{
    return 0;
}

/*
 * \brief Identify differential variables
 */
template <class ScalarT, typename IdxT>
int DiscreteGenerator<ScalarT, IdxT>::tagDifferentiable()
{
    return 0;
}

/**
 * @brief Contributes to the bus residual.
 *
 * Must be connected to a PQ bus.
 */
template <class ScalarT, typename IdxT>
int DiscreteGenerator<ScalarT, IdxT>::evaluateResidual()
{
    //   ### Externals Componenets ###
    //Reference P
    ScalarT wcom = this->wb_ - this->refmp_ * this->y_[0];

    this->f_[0] = 0;
    this->f_[1] = 0;
    this->f_[2] = 0;

    // ### Internal Componenets ##
    f_[3] = -yp_[3] + wb_ - mp_ * y_[4] - wcom;

    f_[4] = -yp_[4] + wc_ * (y_[12] * y_[14] + y_[13] * y_[15] - y_[4]);
    f_[5] =  -yp_[4] + wc_ * (-y_[12] * y_[15] + y_[13] * y_[14] - y_[5]);

    ScalarT vod_star = Vn_ - nq_ * y_[5];
    ScalarT voq_star = 0;

    f_[6] = -yp_[6] + vod_star - y_[12];
    f_[7] = -yp_[7] + voq_star - y_[13];

    ScalarT ild_star = F_ * y_[14] - wb_ * Cf_ * y_[13] + Kpv_ * (vod_star - y_[12]) + Kiv_ * y_[6] ;
    ScalarT ilq_star = F_ * y_[15] + wb_ * Cf_ * y_[12] + Kpv_ * (voq_star - y_[13]) + Kiv_ * y_[7];

    f_[8] = -yp_[8] + ild_star - y_[10];
    f_[9] = -yp_[9] + ilq_star - y_[11];

    ScalarT vid_star = -wb_ * Lf_ * y_[11] + Kpc_ * (ild_star - y_[10]) + Kic_ * y_[8];
    ScalarT viq_star = wb_ * Lf_ * y_[10] + Kpc_ * (ilq_star - y_[11]) + Kic_ * y_[9]; 
    
    f_[10] = -yp_[10] - (rLf_ / Lf_) * y_[10] + wcom * y_[11] + (1/Lf_) * (vid_star - y_[12]);
    f_[11] = -yp_[11] - (rLf_ / Lf_) * y_[11] - wcom * y_[10] + (1/Lf_) * (viq_star - y_[13]);

    f_[12] = -yp_[12] + wcom * y_[13] + (1/Cf_) * (y_[10] - y_[14]);
    f_[13] = -yp_[13] - wcom * y_[12] + (1/Cf_) * (y_[11] - y_[15]);

    f_[14] = -yp_[14] - (rLc_ / Lc_) * y_[14] + wcom * y_[15] + (1/Lc_)*(y_[12] - y_[1]);
    f_[15] = -yp_[15] - (rLc_ / Lc_) * y_[15] - wcom * y_[14] + (1/Lc_)*(y_[13] - y_[2]);
    return 0;
}

/**
 * @brief Compute the jacobian of the DiscreteGenerator for iteration. dF/dy - \alpha dF/dy'
 * 
 * The matrix dF/dy should be
 * [         0,     0,     0, 0,   0,                0,            0,            0,      0,      0,                 0,                 0,                 0,                 0,             0,             0]
[         0,     0,     0, 0,   0,                0,            0,            0,      0,      0,                 0,                 0,                 0,                 0,             0,             0]
[         0,     0,     0, 0,   0,                0,            0,            0,      0,      0,                 0,                 0,                 0,                 0,             0,             0]
[     mpref,     0,     0, 0, -mp,                0,            0,            0,      0,      0,                 0,                 0,                 0,                 0,             0,             0]
[         0,     0,     0, 0, -wc,                0,            0,            0,      0,      0,                 0,                 0,            wc*x15,            wc*x16,        wc*x13,        wc*x14]
[         0,     0,     0, 0,   0,              -wc,            0,            0,      0,      0,                 0,                 0,           -wc*x16,            wc*x15,        wc*x14,       -wc*x13]
[         0,     0,     0, 0,   0,              -nq,            0,            0,      0,      0,                 0,                 0,                -1,                 0,             0,             0]
[         0,     0,     0, 0,   0,                0,            0,            0,      0,      0,                 0,                 0,                 0,                -1,             0,             0]
[         0,     0,     0, 0,   0,          -Kpv*nq,          Kiv,            0,      0,      0,                -1,                 0,              -Kpv,            -Cf*wb,             F,             0]
[         0,     0,     0, 0,   0,                0,            0,          Kiv,      0,      0,                 0,                -1,             Cf*wb,              -Kpv,             0,             F]
[-mpref*x12,     0,     0, 0,   0, -(Kpc*Kpv*nq)/Lf, (Kiv*Kpc)/Lf,            0, Kic/Lf,      0, - Kpc/Lf - rLf/Lf,         -mpref*x1, -(Kpc*Kpv + 1)/Lf,   -(Cf*Kpc*wb)/Lf,    (F*Kpc)/Lf,             0]
[ mpref*x11,     0,     0, 0,   0,                0,            0, (Kiv*Kpc)/Lf,      0, Kic/Lf,          mpref*x1, - Kpc/Lf - rLf/Lf,    (Cf*Kpc*wb)/Lf, -(Kpc*Kpv + 1)/Lf,             0,    (F*Kpc)/Lf]
[-mpref*x14,     0,     0, 0,   0,                0,            0,            0,      0,      0,              1/Cf,                 0,                 0,     wb - mpref*x1,         -1/Cf,             0]
[ mpref*x13,     0,     0, 0,   0,                0,            0,            0,      0,      0,                 0,              1/Cf,     mpref*x1 - wb,                 0,             0,         -1/Cf]
[-mpref*x16, -1/Lc,     0, 0,   0,                0,            0,            0,      0,      0,                 0,                 0,              1/Lc,                 0,       -rLc/Lc, wb - mpref*x1]
[ mpref*x15,     0, -1/Lc, 0,   0,                0,            0,            0,      0,      0,                 0,                 0,                 0,              1/Lc, mpref*x1 - wb,       -rLc/Lc]
 * 'Generated from MATLAB symbolic'
 * Jacobian is mostly constant besides reference and rows 4 & 5
 * 
 * @tparam ScalarT 
 * @tparam IdxT 
 * @return int 
 */
template <class ScalarT, typename IdxT>
int DiscreteGenerator<ScalarT, IdxT>::evaluateJacobian()
{
    //Create dF/dy'
    std::vector<IdxT> rcordder(13);
    std::vector<ScalarT> valsder(13,1.0);
    for (int i = 0; i < 13; i++)
    {
        rcordder[i] = i + 3;
    }
    COO_Matrix<ScalarT,IdxT> Jacder = COO_Matrix<ScalarT, IdxT>(rcordder, rcordder, valsder,16,16);

    //Create dF/dy
    //r = 3
    std::vector<IdxT> ctemp{0, 4};
    std::vector<IdxT> rtemp(ctemp.size(),3);
    std::vector<ScalarT> valtemp{this->refmp_, -this->mp_};
    this->J_.setValues(rtemp, ctemp, valtemp);

    //r = 4
    ctemp = {4, 12, 13, 14, 15};
    std::fill(rtemp.begin(), rtemp.end(), 4);
    valtemp = {-this->wc_, this->wc_*y_[14], this->wc_*y_[15], this->wc_*y_[12], this->wc_*y_[13]};
    this->J_.setValues(rtemp, ctemp, valtemp);

    //r = 5
    ctemp = {5, 12, 13, 14, 15};
    std::fill(rtemp.begin(), rtemp.end(), 5);
    valtemp = {-this->wc_, -this->wc_*y_[15], this->wc_*y_[14], this->wc_*y_[13], -this->wc_*y_[12]};
    this->J_.setValues(rtemp, ctemp, valtemp);

    //r = 6
    ctemp = {5, 12};
    std::fill(rtemp.begin(), rtemp.end(), 6);
    valtemp = {-this->nq_, -1.0};
    this->J_.setValues(rtemp, ctemp, valtemp);
    
    //r = 7
    ctemp = {13};
    std::fill(rtemp.begin(), rtemp.end(), 7);
    valtemp = {-1.0};
    this->J_.setValues(rtemp, ctemp, valtemp);

    //r = 8
    ctemp = {5,6,10,12,13,14};
    std::fill(rtemp.begin(), rtemp.end(), 8);
    valtemp = {-this->Kpv_*this->nq_, this->Kiv_, -1.0, -this->Kpv_, -this->Cf_*this->wb_, this->F_};
    this->J_.setValues(rtemp, ctemp, valtemp);

    //r = 9
    ctemp = {7, 11, 12, 13, 15};
    std::fill(rtemp.begin(), rtemp.end(), 9);
    valtemp = {this->Kiv_, -1.0, this->Cf_*this->wb_,-this->Kpv_,this->F_};
    this->J_.setValues(rtemp, ctemp, valtemp);
    
    //r = 10
    ctemp = {0, 5, 6, 8, 10, 11, 12, 13, 14};
    std::fill(rtemp.begin(), rtemp.end(), 10);
    valtemp = {-this->refmp_ * y_[11], -(this->Kpc_ * this->Kpv_ * this->nq_) / this->Lf_, (this->Kpc_ * this->Kiv_) / this->Lf_, this->Kic_ / this->Lf_, -(this->Kpc_ + this->rLf_) / this->Lf_, -this->refmp_ * y_[0], -(this->Kpc_ * this->Kpv_ + 1.0) / this->Lf_, -(this->Cf_ * this->Kpc_ * this->wb_) / this->Lf_, (this->F_ * this->Kpc_) / this->Lf_};
    this->J_.setValues(rtemp, ctemp, valtemp);

    //r = 11
    ctemp = {0, 7, 9, 10, 11, 12, 13, 15};
    std::fill(rtemp.begin(), rtemp.end(), 11);
    valtemp = {this->refmp_ * y_[10], (this->Kiv_ * this->Kpc_) / this->Lf_, this->Kic_ / this->Lf_, this->refmp_ * y_[0], -(this->Kpc_ + this->rLf_) / this->Lf_, (this->Cf_ * this->Kpc_ * this->wb_) / this->Lf_, -(this->Kpc_ * this->Kpv_ + 1.0) / this->Lf_, (this->F_ * this->Kpc_) / this->Lf_};
    this->J_.setValues(rtemp, ctemp, valtemp);

    //r = 12
    ctemp = {0, 10, 13, 14};
    std::fill(rtemp.begin(), rtemp.end(), 12);
    valtemp = {-this->refmp_ * y_[13], 1.0 / this->Cf_, this->wb_ - this->refmp_ * y_[0], -1.0 / this->Cf_};
    this->J_.setValues(rtemp, ctemp, valtemp);

    
    //r = 13
    ctemp = {0, 11, 12, 15};
    std::fill(rtemp.begin(), rtemp.end(), 13);
    valtemp = {this->refmp_ * y_[12], 1.0 / this->Cf_, -this->wb_ + this->refmp_ * y_[0], -1.0 / this->Cf_};
    this->J_.setValues(rtemp, ctemp, valtemp);

    
    //r = 14
    ctemp = {0, 1, 12, 14, 15};
    std::fill(rtemp.begin(), rtemp.end(), 14);
    valtemp = {-this->refmp_ * y_[15], -1.0 / this->Lc_, 1.0 / this->Lc_, -this->rLc_ / this->Lc_, this->wb_ - this->refmp_ * y_[0]};
    this->J_.setValues(rtemp, ctemp, valtemp);

    
    //r = 15
    ctemp = {0, 2, 13, 14, 15};
    std::fill(rtemp.begin(), rtemp.end(), 15);
    valtemp = {-this->refmp_ * y_[14], -1.0 / this->Lc_, 1.0 / this->Lc_, -this->wb_ + this->refmp_ * y_[0], -this->rLc_ / this->Lc_};
    this->J_.setValues(rtemp, ctemp, valtemp);


    //Perform dF/dy + \alpha dF/dy'

    this->J_.AXPY(-this->alpha_, Jacder);

    return 0;
}

template <class ScalarT, typename IdxT>
int DiscreteGenerator<ScalarT, IdxT>::evaluateIntegrand()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int DiscreteGenerator<ScalarT, IdxT>::initializeAdjoint()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int DiscreteGenerator<ScalarT, IdxT>::evaluateAdjointResidual()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int DiscreteGenerator<ScalarT, IdxT>::evaluateAdjointIntegrand()
{
    return 0;
}




// Available template instantiations
template class DiscreteGenerator<double, long int>;
template class DiscreteGenerator<double, size_t>;


} //namespace ModelLib

