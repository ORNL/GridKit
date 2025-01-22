
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
 * 
 * Model is from paper: "Modeling, Analysis and Testing of Autonomous Operation
 * of an Inverter-Based Microgrid", Nagaraju Pogaku, Milan Prodanovic, and
 * Timothy C. Green, Section C
 * 
 * @todo Consider having \omegaref as a global constant, not a node variable.
 */

template <class ScalarT, typename IdxT>
MicrogridLine<ScalarT, IdxT>::MicrogridLine(IdxT id, ScalarT R,ScalarT L)
  : R_(R),
    L_(L)
{
    // internals [id, iq]
    // externals [\omegaref, vbd_in, vbq_in, vbd_out, vbq_out]
    size_ = 7;
    n_intern_ = 2;
    n_extern_ = 5;
    extern_indices_ = {0,1,2,3,4};
    idc_ = id;
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
    y_.resize(size_);
    yp_.resize(size_);
    f_.resize(size_);
    
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
 * 
 */
template <class ScalarT, typename IdxT>
int MicrogridLine<ScalarT, IdxT>::evaluateResidual()
{
    //ref motor
    f_[0] = 0.0;

    //input
    f_[1] = -y_[5] ;
    f_[2] = -y_[6] ;

    //output
    f_[3] = y_[5] ;
    f_[4] = y_[6] ;

    //Internal variables
    f_[5] = -yp_[5] - (R_ / L_) * y_[5] + y_[0]*y_[6] + (y_[1] - y_[3])/L_;
    f_[6] = -yp_[6] - (R_ / L_) * y_[6] - y_[0]*y_[5] + (y_[2] - y_[4])/L_;


    return 0;
}

/**
 * @brief Generate Jacobian for Microgrid Line
 * 
 * @tparam ScalarT 
 * @tparam IdxT 
 * @return int 
 */
template <class ScalarT, typename IdxT>
int MicrogridLine<ScalarT, IdxT>::evaluateJacobian()
{
    J_.zeroMatrix();

    //Create dF/dy
    std::vector<IdxT> rtemp{1,2,3,4};
    std::vector<IdxT> ctemp{5,6,5,6};
    std::vector<ScalarT> vals{-1.0,-1.0,1.0,1.0};
    J_.setValues(rtemp, ctemp, vals);

    
    std::vector<IdxT> ccord{0, 1, 3, 5, 6};

    std::vector<IdxT> rcord(ccord.size(),5);
    vals = {y_[6], (1.0 / L_) , -(1.0 / L_),  - (R_ / L_) , y_[0]};
    J_.setValues(rcord, ccord, vals);

    
    std::vector<IdxT> ccor2{0, 2, 4, 5, 6};
    std::fill(rcord.begin(), rcord.end(), 6);
    vals = {-y_[5], (1.0 / L_) , -(1.0 / L_), -y_[0], - (R_ / L_)};
    J_.setValues(rcord, ccor2, vals);


    //Create -dF/dy'
    std::vector<IdxT> rcordder{5,6};
    std::vector<IdxT> ccordder{5,6};
    std::vector<ScalarT> valsder{-1.0, -1.0};
    COO_Matrix<ScalarT,IdxT> Jacder = COO_Matrix<ScalarT, IdxT>(rcordder, ccordder, valsder,7,7);
    
    //Perform dF/dy + \alpha dF/dy'
    J_.axpy(alpha_, Jacder);


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

