


#include <iostream>
#include <cmath>
#include <vector>
#include "LinearTransformer.hpp"

namespace GridKit {

/*!
 * @brief Constructor for a LinearTransformer model
 *
 * Calls default ModelEvaluatorImpl constructor.
 * @todo Not tested in any model yet. Should be
 * @todo Has not been tested for correctness
 * 
 * @tparam ScalarT - floating point type for the model
 * @tparam IdxT - integer index type for the model
 * 
 * @param[in] id - unique identifier for the component
 * @param[in] L0 - inductance 0
 * @param[in] L1 - inductance 1
 * @param[in] R0 - resistance 0
 * @param[in] R1 - resistance 1
 * @param[in] M - mutual inductance
 */

template <class ScalarT, typename IdxT>
LinearTransformer<ScalarT, IdxT>::LinearTransformer(IdxT id, ScalarT L0, ScalarT L1, ScalarT R0, ScalarT R1, ScalarT M)
  : L0_(L0),
    L1_(L1),
    R0_(R0),
    R1_(R1),
    M_(M)
{
    size_ = 4;
    n_intern_ = 2;
    n_extern_ = 2;
    extern_indices_ = {0,1};
    idc_ = id;
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
    y_.resize(size_);
    yp_.resize(size_);
    f_.resize(size_);
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
 * @brief Computes the component resisdual
 * 
 */
template <class ScalarT, typename IdxT>
int LinearTransformer<ScalarT, IdxT>::evaluateResidual()
{
    f_[0] = y_[2];
    f_[1] = y_[3];
    f_[2] = y_[0] - R0_ * y_[2] - L0_ * yp_[2] - M_ * yp_[3];
    f_[2] = y_[1] - R1_ * y_[3] - M_ * yp_[2] - L1_ * yp_[3];
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


} //namespace GridKit
