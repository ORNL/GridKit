


#include <iostream>
#include <cmath>
#include <vector>
#include "VoltageSource.hpp"

namespace ModelLib {

/*!
 * @brief Constructor for a constant VoltageSource model
 *
 * Calls default ModelEvaluatorImpl constructor.
 */

template <class ScalarT, typename IdxT>
VoltageSource<ScalarT, IdxT>::VoltageSource(IdxT id, ScalarT V)
  : V_(V)
{
	this->size_ = 3;
    this->n_intern = 1;
    this->n_extern = 2;
    this->extern_indices = {0,1};
    this->idc_ = id;
}

template <class ScalarT, typename IdxT>
VoltageSource<ScalarT, IdxT>::~VoltageSource()
{
}

/*!
 * @brief allocate method computes sparsity pattern of the Jacobian.
 */
template <class ScalarT, typename IdxT>
int VoltageSource<ScalarT, IdxT>::allocate()
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
int VoltageSource<ScalarT, IdxT>::initialize()
{
    return 0;
}

/*
 * \brief Identify differential variables
 */
template <class ScalarT, typename IdxT>
int VoltageSource<ScalarT, IdxT>::tagDifferentiable()
{
    return 0;
}

/**
 * @brief Contributes to the bus residual.
 *
 * Must be connected to a PQ bus.
 */
template <class ScalarT, typename IdxT>
int VoltageSource<ScalarT, IdxT>::evaluateResidual()
{
	//Note this leaves induction lumped into y. Perhaps would be better to seperate volatge and induction into seperate vectors
	// for easier development
    this->f_[0] = this->y_[2];
    this->f_[1] = -this->y_[2];
	this->f_[2] = this->y_[1] - this->y_[0] - this->V_;
    return 0;
}

template <class ScalarT, typename IdxT>
int VoltageSource<ScalarT, IdxT>::evaluateJacobian()
{
    
    //Create dF/dy
    std::vector<IdxT> rcord{0,1,2,2};
    std::vector<IdxT> ccord{2,2,0,1};
    std::vector<ScalarT> vals{1.0, -1.0, -1.0, 1.0};
    this->J_.setValues(rcord, ccord, vals);

    return 0;
}

template <class ScalarT, typename IdxT>
int VoltageSource<ScalarT, IdxT>::evaluateIntegrand()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int VoltageSource<ScalarT, IdxT>::initializeAdjoint()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int VoltageSource<ScalarT, IdxT>::evaluateAdjointResidual()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int VoltageSource<ScalarT, IdxT>::evaluateAdjointIntegrand()
{
    return 0;
}





// Available template instantiations
template class VoltageSource<double, long int>;
template class VoltageSource<double, size_t>;


} //namespace ModelLib

