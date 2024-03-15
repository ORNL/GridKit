


#include <iostream>
#include <cmath>
#include <vector>
#include "SynchronousMachine.hpp"


namespace ModelLib {



/*!
 * @brief Constructor for a constant SynchronousMachine model
 *
 * Calls default ModelEvaluatorImpl constructor.
 */

template <class ScalarT, typename IdxT>
SynchronousMachine<ScalarT, IdxT>::SynchronousMachine(IdxT id, ScalarT Lls, std::tuple<ScalarT,ScalarT> Llkq, ScalarT Llfd, ScalarT Llkd, ScalarT Lmq, ScalarT Lmd, ScalarT Rs, std::tuple<ScalarT,ScalarT> Rkq, ScalarT Rfd, ScalarT Rkd, ScalarT J, ScalarT P, ScalarT mub)
  : Lls_(Lls),
	Llkq_(Llkq),
	Llfd_(Llfd),
	Llkd_(Llkd),
	Lmq_(Lmq),
	Lmd_(Lmd),
	Rs_(Rs),
	Rkq_(Rkq),
	Rfd_(Rfd),
	Rkd_(Rkd),
	J_(J),
	P_(P),
	mub_(mub)
{
	this->size_ = 13;
    this->n_intern = 6;
    this->n_extern = 7;
    this->extern_indices = {0,1,2,3,4};
    this->idc_ = id;
}

template <class ScalarT, typename IdxT>
SynchronousMachine<ScalarT, IdxT>::~SynchronousMachine()
{
}

/*!
 * @brief allocate method computes sparsity pattern of the Jacobian.
 */
template <class ScalarT, typename IdxT>
int SynchronousMachine<ScalarT, IdxT>::allocate()
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
int SynchronousMachine<ScalarT, IdxT>::initialize()
{
    return 0;
}

/*
 * \brief Identify differential variables
 */
template <class ScalarT, typename IdxT>
int SynchronousMachine<ScalarT, IdxT>::tagDifferentiable()
{
    return 0;
}

/**
 * @brief Contributes to the bus residual.
 *
 * Must be connected to a PQ bus.
 */
template <class ScalarT, typename IdxT>
int SynchronousMachine<ScalarT, IdxT>::evaluateResidual()
{
	ScalarT rkq1 = std::get<0>(Rkq_);
	ScalarT rkq2 = std::get<1>(Rkq_);
	ScalarT llkq1 = std::get<0>(Llkq_);
	ScalarT llkq2 = std::get<1>(Llkq_);

	ScalarT cos1 = cos((P_/2.0)*y_[5]);
	ScalarT sin1 = sin((P_/2.0)*y_[5]);
	ScalarT cos23m = cos((P_/2.0)*y_[5] - (2.0/3.0)*M_PI);
	ScalarT sin23m = sin((P_/2.0)*y_[5] - (2.0/3.0)*M_PI);
	ScalarT cos23p = cos((P_/2.0)*y_[5] + (2.0/3.0)*M_PI);
	ScalarT sin23p = sin((P_/2.0)*y_[5] + (2.0/3.0)*M_PI);

    this->f_[0] = y_[6]*cos1 + y_[7]*sin1 + y_[8];
    this->f_[1] = y_[6]*cos23m + y_[7]*sin23m + y_[8];
	this->f_[2] = y_[6]*cos23p + y_[7]*sin23p + y_[8];
	this->f_[3] = J_ * yp_[4] - (3.0/4.0)*P_*(Lmd_ *y_[6]* (y_[7] + y_[11] + y_[12]) - Lmq_*y_[7]*(y_[6] + y_[9] + y_[0]));
	this->f_[4] = yp_[5] - y_[4];
	this->f_[5] = (-2.0/3.0)*(y_[0]*cos1 +y_[1]*cos23m + y_[2]*cos23p) + Rs_*y_[6] + (Lls_ + Lmq_)*yp_[6] + Lmq_*yp_[9] + Lmq_*yp_[10] + y_[4]*(P_/2.0)*((Lls_ + Lmd_)*y_[7] + Lmd_*y_[11] + Lmd_*y_[12]);
	this->f_[6] = (-2.0/3.0)*(y_[0]*sin1 -y_[1]*sin23m - y_[2]*sin23p) + Rs_*y_[7] + (Lls_ + Lmd_)*yp_[7] + Lmd_*yp_[11] + Lmd_*yp_[12] - y_[4]*(P_/2.0)*((Lls_ + Lmq_)*y_[6] + Lmq_*y_[9] + Lmq_*y_[10]);
	this->f_[7] = (-1.0/3.0)*(y_[0] + y_[1] + y_[2]) + Rs_*y_[8] + Lls_*yp_[8];
	this->f_[8] = rkq1*y_[9] + (llkq1 + Lmq_)*yp_[9] + Lmq_*yp_[6] + Lmq_*yp_[10];
	this->f_[9] = rkq1*y_[9] + (llkq1 + Lmq_)*yp_[9] + Lmq_*yp_[6] + Lmq_*yp_[10];
    return 0;
}

template <class ScalarT, typename IdxT>
int SynchronousMachine<ScalarT, IdxT>::evaluateJacobian()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int SynchronousMachine<ScalarT, IdxT>::evaluateIntegrand()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int SynchronousMachine<ScalarT, IdxT>::initializeAdjoint()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int SynchronousMachine<ScalarT, IdxT>::evaluateAdjointResidual()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int SynchronousMachine<ScalarT, IdxT>::evaluateAdjointIntegrand()
{
    return 0;
}





// Available template instantiations
template class SynchronousMachine<double, long int>;
template class SynchronousMachine<double, size_t>;


} //namespace ModelLib

