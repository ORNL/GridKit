


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
    //input
    this->f_[0] = C_ * this->yp_[2];
    //output
    this->f_[1] = -C_ * this->yp_[2];

    //internal
    this->f_[2] =  -C_ * this->yp_[2] + this->y_[0] - this->y_[1] - this->y_[2];
    return 0;
}

template <class ScalarT, typename IdxT>
int Capacitor<ScalarT, IdxT>::evaluateJacobian()
{
    this->J_.zeroMatrix();
    //Create dF/dy
    std::vector<IdxT> rcord{2,2,2};
    std::vector<IdxT> ccord{0,1,2};
    std::vector<ScalarT> vals{1.0, -1.0, -1.0};
    this->J_.setValues(rcord, ccord, vals);

    //Create -dF/dy'
    std::vector<IdxT> rcordder{0,1,2};
    std::vector<IdxT> ccordder{2,2,2};
    std::vector<ScalarT> valsder{C_, -C_, -this->C_};
    COO_Matrix<ScalarT,IdxT> Jacder = COO_Matrix<ScalarT, IdxT>(rcordder, ccordder, valsder,3,3);
    
    //Perform dF/dy + \alpha dF/dy'
    this->J_.AXPY(this->alpha_, Jacder);

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
