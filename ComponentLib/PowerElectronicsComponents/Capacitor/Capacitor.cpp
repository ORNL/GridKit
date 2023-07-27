


#include <iostream>
#include <cmath>
#include <vector>
#include "Capacitor.hpp"

namespace ModelLib {

/*!
 * @brief Constructor for a constant load model
 *
 * Calls default ModelEvaluatorImpl constructor.
 */

template <class ScalarT, typename IdxT>
Capacitor<ScalarT, IdxT>::Capacitor(IdxT id, ScalarT C)
  :  C_(C)
{
	this->size_ = 3;
    this->n_intern = 1;
    this->n_extern = 2;
    this->extern_indices = {0,1};
    this->idc_ = id;
}

template <class ScalarT, typename IdxT>
Capacitor<ScalarT, IdxT>::~Capacitor()
{
}

/*!
 * @brief allocate method computes sparsity pattern of the Jacobian.
 */
template <class ScalarT, typename IdxT>
int Capacitor<ScalarT, IdxT>::allocate()
{
    this->y_.resize(this->size_);
    this->yp_.resize(this->size_);
    this->f_.resize(this->size_);
    this->J_.resize(this->size_^2);
    this->M_.resize(this->size_^2);
    return 0;
}

/**
 * Initialization of the grid model
 */
template <class ScalarT, typename IdxT>
int Capacitor<ScalarT, IdxT>::initialize()
{
    return 0;
}

/*
 * \brief Identify differential variables
 */
template <class ScalarT, typename IdxT>
int Capacitor<ScalarT, IdxT>::tagDifferentiable()
{
    return 0;
}

/**
 * @brief Contributes to the bus residual.
 *
 * Must be connected to a PQ bus.
 */
template <class ScalarT, typename IdxT>
int Capacitor<ScalarT, IdxT>::evaluateResidual()
{
    this->f_[0] = this->yp_[2];
    this->f_[1] = -this->yp_[2];
    this->f_[2] = this->y_[0] - this->y_[1] - this->y_[2] - this->C_ * this->yp_[2];
    return 0;
}

template <class ScalarT, typename IdxT>
int Capacitor<ScalarT, IdxT>::evaluateJacobian()
{

    return 0;
}

template <class ScalarT, typename IdxT>
int Capacitor<ScalarT, IdxT>::evaluateIntegrand()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int Capacitor<ScalarT, IdxT>::initializeAdjoint()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int Capacitor<ScalarT, IdxT>::evaluateAdjointResidual()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int Capacitor<ScalarT, IdxT>::evaluateAdjointIntegrand()
{
    return 0;
}




// Available template instantiations
template class Capacitor<double, long int>;
template class Capacitor<double, size_t>;


} //namespace ModelLib

