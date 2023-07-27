


#include <iostream>
#include <cmath>
#include <vector>
#include "LinearTransformer.hpp"

namespace ModelLib {

/*!
 * @brief Constructor for a constant LinearTransformer model
 *
 * Calls default ModelEvaluatorImpl constructor.
 */

template <class ScalarT, typename IdxT>
LinearTransformer<ScalarT, IdxT>::LinearTransformer(IdxT id, ScalarT L1, ScalarT L2, ScalarT R1, ScalarT R2, ScalarT M)
  : L1_(L1),
	L2_(L2),
	R1_(R1),
	R2_(R2),
	M_(M)
{
	this->size_ = 4;
    this->n_intern = 2;
    this->n_extern = 2;
    this->extern_indices = {0,1};
    this->idc_ = id;
}

template <class ScalarT, typename IdxT>
LinearTransformer<ScalarT, IdxT>::~LinearTransformer()
{
}

/*!
 * @brief allocate method computes sparsity pattern of the Jacobian.
 */
template <class ScalarT, typename IdxT>
int LinearTransformer<ScalarT, IdxT>::allocate()
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
int LinearTransformer<ScalarT, IdxT>::initialize()
{
    return 0;
}

/*
 * \brief Identify differential variables
 */
template <class ScalarT, typename IdxT>
int LinearTransformer<ScalarT, IdxT>::tagDifferentiable()
{
    return 0;
}

/**
 * @brief Contributes to the bus residual.
 *
 * Must be connected to a PQ bus.
 */
template <class ScalarT, typename IdxT>
int LinearTransformer<ScalarT, IdxT>::evaluateResidual()
{
	//Note this leaves induction lumped into y. Perhaps would be better to seperate volatge and induction into seperate vectors
	// for easier development
    this->f_[0] = this->y_[2];
    this->f_[1] = this->y_[3];
	this->f_[2] = this->y_[0] - this->R1_ * this->y_[2] - this->L1_ * this->yp_[2] - this->M_ * this->yp_[3];
	this->f_[2] = this->y_[1] - this->R2_ * this->y_[3] - this->M_ * this->yp_[2] - this->L2_ * this->yp_[3];
    return 0;
}

template <class ScalarT, typename IdxT>
int LinearTransformer<ScalarT, IdxT>::evaluateJacobian()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int LinearTransformer<ScalarT, IdxT>::evaluateIntegrand()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int LinearTransformer<ScalarT, IdxT>::initializeAdjoint()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int LinearTransformer<ScalarT, IdxT>::evaluateAdjointResidual()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int LinearTransformer<ScalarT, IdxT>::evaluateAdjointIntegrand()
{
    return 0;
}





// Available template instantiations
template class LinearTransformer<double, long int>;
template class LinearTransformer<double, size_t>;


} //namespace ModelLib

