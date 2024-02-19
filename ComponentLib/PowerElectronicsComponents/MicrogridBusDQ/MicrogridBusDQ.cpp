
#include <iostream>
#include <cmath>
#include <vector>
#include "MicrogridBusDQ.hpp"

namespace ModelLib {

/*!
 * @brief Constructor for a constant MicrogridBusDQ model
 *
 * Calls default ModelEvaluatorImpl constructor.
 * 
 * Each microgrid line has a virtual resistance RN
 */

template <class ScalarT, typename IdxT>
MicrogridBusDQ<ScalarT, IdxT>::MicrogridBusDQ(IdxT id, ScalarT RN)
  : RN_(RN) 
{
    // externals [vbus_d, vbus_q]
	this->size_ = 2;
    this->n_intern = 0;
    this->n_extern = 2;
    this->extern_indices = {0,1};
    this->idc_ = id;
}

template <class ScalarT, typename IdxT>
MicrogridBusDQ<ScalarT, IdxT>::~MicrogridBusDQ()
{
}

/*!
 * @brief allocate method computes sparsity pattern of the Jacobian.
 */
template <class ScalarT, typename IdxT>
int MicrogridBusDQ<ScalarT, IdxT>::allocate()
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
int MicrogridBusDQ<ScalarT, IdxT>::initialize()
{
    return 0;
}

/*
 * \brief Identify differential variables
 */
template <class ScalarT, typename IdxT>
int MicrogridBusDQ<ScalarT, IdxT>::tagDifferentiable()
{
    return 0;
}

/**
 * @brief Evaluate residual of microgrid line
 * This model has "Virtual resistors" on both ends with parameters RN1 and RN2
 *
 */
template <class ScalarT, typename IdxT>
int MicrogridBusDQ<ScalarT, IdxT>::evaluateResidual()
{
    //bus voltage
    f_[0] = -y_[0] / RN_;
    f_[1] = -y_[1] / RN_;

    return 0;
}

/**
 * @brief Generate Jacobian for Transmission Line
 * 
 * @tparam ScalarT 
 * @tparam IdxT 
 * @return int 
 */
template <class ScalarT, typename IdxT>
int MicrogridBusDQ<ScalarT, IdxT>::evaluateJacobian()
{
    this->J_.zeroMatrix();

    //Create dF/dy
    std::vector<IdxT> rtemp{0,1};
    std::vector<IdxT> ctemp{0,1};
    std::vector<ScalarT> vals{-1.0 / RN_,-1.0 / RN_};
    this->J_.setValues(rtemp, ctemp, vals);

    return 0;
}

template <class ScalarT, typename IdxT>
int MicrogridBusDQ<ScalarT, IdxT>::evaluateIntegrand()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int MicrogridBusDQ<ScalarT, IdxT>::initializeAdjoint()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int MicrogridBusDQ<ScalarT, IdxT>::evaluateAdjointResidual()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int MicrogridBusDQ<ScalarT, IdxT>::evaluateAdjointIntegrand()
{
    return 0;
}


// Available template instantiations
template class MicrogridBusDQ<double, long int>;
template class MicrogridBusDQ<double, size_t>;


} //namespace ModelLib

