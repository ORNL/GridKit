


#include <iostream>
#include <cmath>
#include <vector>
#include "Capacitor.hpp"

namespace ModelLib {

/*!
 * @brief Constructor for Capacitor
 *
 * @todo this needs to be tested on some circuit
 *
 * Calls default ModelEvaluatorImpl constructor.
 */

template <class ScalarT, typename IdxT>
Capacitor<ScalarT, IdxT>::Capacitor(IdxT id, ScalarT C)
  :  C_(C)
{
    size_ = 3;
    n_intern_ = 1;
    n_extern_ = 2;
    extern_indices_ = {0,1};
    idc_ = id;
}

template <class ScalarT, typename IdxT>
Capacitor<ScalarT, IdxT>::~Capacitor()
{
}

/*!
 * @brief allocate method creates memory for vectors
 */
template <class ScalarT, typename IdxT>
int Capacitor<ScalarT, IdxT>::allocate()
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
 * @brief Evaluate the resisdual of the Capcitor
 *
 */
template <class ScalarT, typename IdxT>
int Capacitor<ScalarT, IdxT>::evaluateResidual()
{
    //input
    f_[0] = C_ * yp_[2];
    //output
    f_[1] = -C_ * yp_[2];

    //internal
    f_[2] =  -C_ * yp_[2] + y_[0] - y_[1] - y_[2];
    return 0;
}

/**
 * @brief Compute the Jacobian dF/dy - a dF/dy'
 * 
 * @tparam ScalarT 
 * @tparam IdxT 
 * @return int 
 */
template <class ScalarT, typename IdxT>
int Capacitor<ScalarT, IdxT>::evaluateJacobian()
{
    J_.zeroMatrix();
    //Create dF/dy
    std::vector<IdxT> rcord{2,2,2};
    std::vector<IdxT> ccord{0,1,2};
    std::vector<ScalarT> vals{1.0, -1.0, -1.0};
    J_.setValues(rcord, ccord, vals);

    //Create dF/dy'
    std::vector<IdxT> rcordder{0,1,2};
    std::vector<IdxT> ccordder{2,2,2};
    std::vector<ScalarT> valsder{C_, -C_, -C_};
    COO_Matrix<ScalarT,IdxT> Jacder = COO_Matrix<ScalarT, IdxT>(rcordder, ccordder, valsder,3,3);
    
    //Perform dF/dy + \alpha dF/dy'
    J_.axpy(alpha_, Jacder);

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

