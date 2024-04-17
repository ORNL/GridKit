


#include <iostream>
#include <cmath>
#include <vector>
#include "Resistor.hpp"

namespace ModelLib {

/*!
 * @brief Constructor for a resistor model
 *
 * Calls default ModelEvaluatorImpl constructor.
 */

template <class ScalarT, typename IdxT>
Resistor<ScalarT, IdxT>::Resistor(IdxT id, ScalarT R)
  : R_(R)
{
    size_ = 2;
    n_intern_ = 0;
    n_extern_ = 2;
    extern_indices_ = {0,1};
    idc_ = id;
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
    y_.resize(size_);
    yp_.resize(size_);
    f_.resize(size_);
    
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
 * @brief Computes the resistors resisdual
 *
 */
template <class ScalarT, typename IdxT>
int Resistor<ScalarT, IdxT>::evaluateResidual()
{
    //input
    f_[0] = (y_[0] - y_[1])/R_ ;
    //ouput
    f_[1] = (y_[1] - y_[0])/R_ ;
    return 0;
}

template <class ScalarT, typename IdxT>
int Resistor<ScalarT, IdxT>::evaluateJacobian()
{
    
    //Create dF/dy
    //does compiler make constant???
    std::vector<IdxT> rcord{0,0,1,1};
    std::vector<IdxT> ccord{0,1,0,1};
    std::vector<ScalarT> vals{1.0 / R_, -1.0 / R_, -1.0 / R_, 1.0 / R_};
    J_.setValues(rcord, ccord, vals);

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

