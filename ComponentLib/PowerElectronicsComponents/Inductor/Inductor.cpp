


#include <iostream>
#include <cmath>
#include <vector>
#include "Inductor.hpp"

namespace ModelLib {

/*!
 * @brief Constructor for a inductor
 *
 * Calls default ModelEvaluatorImpl constructor.
 */

template <class ScalarT, typename IdxT>
Inductor<ScalarT, IdxT>::Inductor(IdxT id, ScalarT L)
  : L_(L)
{
    size_ = 3;
    n_intern_ = 1;
    n_extern_ = 2;
    extern_indices_ = {0,1};
    idc_ = id;
}

template <class ScalarT, typename IdxT>
Inductor<ScalarT, IdxT>::~Inductor()
{
}

/*!
 * @brief allocate method computes sparsity pattern of the Jacobian.
 */
template <class ScalarT, typename IdxT>
int Inductor<ScalarT, IdxT>::allocate()
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
int Inductor<ScalarT, IdxT>::initialize()
{
    return 0;
}

/*
 * \brief Identify differential variables
 */
template <class ScalarT, typename IdxT>
int Inductor<ScalarT, IdxT>::tagDifferentiable()
{
    return 0;
}

/**
 * @brief Compute the resisdual of the component
 *
 */
template <class ScalarT, typename IdxT>
int Inductor<ScalarT, IdxT>::evaluateResidual()
{
    //input
    f_[0] = -y_[2];
    //output
    f_[1] = y_[2];
    //internal
    f_[2] = -L_ * yp_[2] + y_[1] - y_[0] ;
    return 0;
}

/**
 * @brief Evaluate the jacobian of the component
 * 
 * @tparam ScalarT 
 * @tparam IdxT 
 * @return int 
 */
template <class ScalarT, typename IdxT>
int Inductor<ScalarT, IdxT>::evaluateJacobian()
{
    J_.zeroMatrix();

    //Create dF/dy
    std::vector<IdxT> rcord{0,1,2,2};
    std::vector<IdxT> ccord{2,2,0,1};
    std::vector<ScalarT> vals{-1.0, 1.0, -1.0, 1.0};
    J_.setValues(rcord, ccord, vals);

    //Create dF/dy'
    std::vector<IdxT> rcordder{2};
    std::vector<IdxT> ccordder{2};
    std::vector<ScalarT> valsder{-L_};
    COO_Matrix<ScalarT,IdxT> Jacder = COO_Matrix<ScalarT, IdxT>(rcordder, ccordder, valsder,3,3);
    
    //Perform dF/dy + \alpha dF/dy'
    J_.axpy(alpha_, Jacder);

    return 0;
}

template <class ScalarT, typename IdxT>
int Inductor<ScalarT, IdxT>::evaluateIntegrand()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int Inductor<ScalarT, IdxT>::initializeAdjoint()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int Inductor<ScalarT, IdxT>::evaluateAdjointResidual()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int Inductor<ScalarT, IdxT>::evaluateAdjointIntegrand()
{
    return 0;
}




// Available template instantiations
template class Inductor<double, long int>;
template class Inductor<double, size_t>;


} //namespace ModelLib

