
#include <iostream>
#include <cmath>
#include <vector>
#include "MicrogridLine.hpp"

namespace ModelLib {

/*!
 * @brief Constructor for a constant MicrogridLine model
 *
 * Calls default ModelEvaluatorImpl constructor.
 * 
 * Each microgrid line has a virtual resistance RN
 */

template <class ScalarT, typename IdxT>
MicrogridLine<ScalarT, IdxT>::MicrogridLine(IdxT id, ScalarT R,ScalarT L)
  : R_(R), L_(L)
{
    // internals [id, iq]
    // externals [\omegaref, vbd_in, vbq_in, vbd_out, vbq_out]
	this->size_ = 7;
    this->n_intern = 2;
    this->n_extern = 5;
    this->extern_indices = {0,1,2,3,4};
    this->idc_ = id;
}

template <class ScalarT, typename IdxT>
MicrogridLine<ScalarT, IdxT>::~MicrogridLine()
{
}

/*!
 * @brief allocate method computes sparsity pattern of the Jacobian.
 */
template <class ScalarT, typename IdxT>
int MicrogridLine<ScalarT, IdxT>::allocate()
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
int MicrogridLine<ScalarT, IdxT>::initialize()
{
    return 0;
}

/*
 * \brief Identify differential variables
 */
template <class ScalarT, typename IdxT>
int MicrogridLine<ScalarT, IdxT>::tagDifferentiable()
{
    return 0;
}

/**
 * @brief Evaluate residual of microgrid line
 * This model has "Virtual resistors" on both ends with parameters RN1 and RN2
 *
 */
template <class ScalarT, typename IdxT>
int MicrogridLine<ScalarT, IdxT>::evaluateResidual()
{
    //ref motor
    this->f_[0] = 0.0;

	//input
    this->f_[1] = -y_[5] ;
    this->f_[2] = -y_[6] ;

    //output
    this->f_[3] = y_[5] ;
    this->f_[4] = y_[6] ;

    //Internal variables
    this->f_[5] = -yp_[5] - (R_ / L_) * y_[5] + y_[0]*y_[6] +  (y_[1] - y_[3])/L_;
    this->f_[6] = -yp_[6] - (R_ / L_) * y_[6] - y_[0]*y_[5] + (y_[2] - y_[4])/L_;


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
int MicrogridLine<ScalarT, IdxT>::evaluateJacobian()
{
    this->J_.zeroMatrix();

    //Create dF/dy
    std::vector<IdxT> rtemp{1,2,3,4};
    std::vector<IdxT> ctemp{5,6,5,6};
    std::vector<ScalarT> vals{-1.0,-1.0,1.0,1.0};
    this->J_.setValues(rtemp, ctemp, vals);

    
    std::vector<IdxT> ccord{0, 1, 3, 5, 6};

    std::vector<IdxT> rcord(ccord.size(),5);
    vals = {y_[6], (1.0 / L_) , -(1.0 / L_),  - (R_ / L_) , y_[0]};
    this->J_.setValues(rcord, ccord, vals);

    
    std::vector<IdxT> ccor2{0, 2, 4, 5, 6};
    std::fill(rcord.begin(), rcord.end(), 6);
    vals = {-y_[5], (1.0 / L_) , -(1.0 / L_), -y_[0], - (R_ / L_)};
    this->J_.setValues(rcord, ccor2, vals);


    //Create -dF/dy'
    std::vector<IdxT> rcordder{5,6};
    std::vector<IdxT> ccordder{5,6};
    std::vector<ScalarT> valsder{-1.0, -1.0};
    COO_Matrix<ScalarT,IdxT> Jacder = COO_Matrix<ScalarT, IdxT>(rcordder, ccordder, valsder,7,7);
    
    //Perform dF/dy + \alpha dF/dy'
    this->J_.AXPY(this->alpha_, Jacder);


    return 0;
}

template <class ScalarT, typename IdxT>
int MicrogridLine<ScalarT, IdxT>::evaluateIntegrand()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int MicrogridLine<ScalarT, IdxT>::initializeAdjoint()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int MicrogridLine<ScalarT, IdxT>::evaluateAdjointResidual()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int MicrogridLine<ScalarT, IdxT>::evaluateAdjointIntegrand()
{
    return 0;
}





// Available template instantiations
template class MicrogridLine<double, long int>;
template class MicrogridLine<double, size_t>;


} //namespace ModelLib

