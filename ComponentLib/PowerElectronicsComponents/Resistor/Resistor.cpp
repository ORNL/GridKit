


#include <iostream>
#include <cmath>
#include <vector>
#include "Resistor.hpp"

namespace ModelLib {

/*!
 * @brief Constructor for a constant resistor model
 *
 * Calls default ModelEvaluatorImpl constructor.
 */

template <class ScalarT, typename IdxT>
Resistor<ScalarT, IdxT>::Resistor(IdxT id, ScalarT R)
  : R_(R)
{
	this->size_ = 2;
    this->n_intern = 0;
    this->n_extern = 2;
    this->extern_indices = {0,1};
    this->idc_ = id;
}

template <class ScalarT, typename IdxT>
Resistor<ScalarT, IdxT>::~Resistor()
{
}

/*!
 * @brief allocate method computes sparsity pattern of the Jacobian.
 */
template <class ScalarT, typename IdxT>
int Resistor<ScalarT, IdxT>::allocate()
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
int Resistor<ScalarT, IdxT>::initialize()
{
    return 0;
}

/*
 * \brief Identify differential variables
 */
template <class ScalarT, typename IdxT>
int Resistor<ScalarT, IdxT>::tagDifferentiable()
{
    return 0;
}

/**
 * @brief Contributes to the bus residual.
 *
 * Must be connected to a PQ bus.
 */
template <class ScalarT, typename IdxT>
int Resistor<ScalarT, IdxT>::evaluateResidual()
{
	
    this->f_[0] = (this->y_[1] - this->y_[0])/this->R_ ;
    this->f_[1] = (this->y_[0] - this->y_[1])/this->R_ ;
    return 0;
}

template <class ScalarT, typename IdxT>
int Resistor<ScalarT, IdxT>::evaluateJacobian()
{
    
    //Create dF/dy
    //does compiler make constant???
    std::vector<IdxT> rcord{0,0,1,1};
    std::vector<IdxT> ccord{0,1,0,1};
    std::vector<ScalarT> vals{-1.0 / this->R_, 1.0 / this->R_, 1.0 / this->R_, -1.0 / this->R_};
    this->J_.setValues(rcord, ccord, vals);

    return 0;
}

template <class ScalarT, typename IdxT>
int Resistor<ScalarT, IdxT>::evaluateIntegrand()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int Resistor<ScalarT, IdxT>::initializeAdjoint()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int Resistor<ScalarT, IdxT>::evaluateAdjointResidual()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int Resistor<ScalarT, IdxT>::evaluateAdjointIntegrand()
{
    return 0;
}





// Available template instantiations
template class Resistor<double, long int>;
template class Resistor<double, size_t>;


} //namespace ModelLib

