
#include <iostream>
#include <cmath>
#include <vector>
#include "MicrogridLoad.hpp"

namespace ModelLib {

/*!
 * @brief Constructor for a constant MicrogridLoad model
 *
 * Calls default ModelEvaluatorImpl constructor.
 * 
 * 
 * Model is from paper: "
    "Modeling, Analysis and Testing of Autonomous Operation of an Inverter-Based Microgrid" Nagaraju Pogaku, Milan Prodanovic, and Timothy C. Green"
 * Section D
 */

template <class ScalarT, typename IdxT>
MicrogridLoad<ScalarT, IdxT>::MicrogridLoad(IdxT id, ScalarT R,ScalarT L)
  : R_(R), 
    L_(L)
{
    // internals [id, iq]
    // externals [\omegaref, vbd_out, vbq_out]
    size_ = 5;
    n_intern_ = 2;
    n_extern_ = 3;
    extern_indices_ = {0,1,2};
    idc_ = id;

}

template <class ScalarT, typename IdxT>
MicrogridLoad<ScalarT, IdxT>::~MicrogridLoad()
{
}

/*!
 * @brief allocate method computes sparsity pattern of the Jacobian.
 */
template <class ScalarT, typename IdxT>
int MicrogridLoad<ScalarT, IdxT>::allocate()
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
int MicrogridLoad<ScalarT, IdxT>::initialize()
{
    return 0;
}

/*
 * \brief Identify differential variables
 */
template <class ScalarT, typename IdxT>
int MicrogridLoad<ScalarT, IdxT>::tagDifferentiable()
{
    return 0;
}

/**
 * @brief Eval Micro Load
 */
template <class ScalarT, typename IdxT>
int MicrogridLoad<ScalarT, IdxT>::evaluateResidual()
{
    //ref motor
    f_[0] = 0.0;

    //only input for loads

    //input
    f_[1] = -y_[3] ;
    f_[2] = -y_[4] ;

    //Internal variables
    f_[3] = -yp_[3] - (R_ / L_) * y_[3] + y_[0]*y_[4] + y_[1] / L_;
    f_[4] = -yp_[4] - (R_ / L_) * y_[4] - y_[0]*y_[3] + y_[2] / L_;


    return 0;
}

/**
 * @brief Generate Jacobian for Micro Load
 * 
 * @tparam ScalarT 
 * @tparam IdxT 
 * @return int 
 */
template <class ScalarT, typename IdxT>
int MicrogridLoad<ScalarT, IdxT>::evaluateJacobian()
{
    J_.zeroMatrix();

    //Create dF/dy
    std::vector<IdxT> rtemp{1,2};
    std::vector<IdxT> ctemp{3,4};
    std::vector<ScalarT> vals{-1.0,-1.0};
    J_.setValues(rtemp, ctemp, vals);

    
    std::vector<IdxT> ccord{0, 1, 3, 4};

    std::vector<IdxT> rcord(ccord.size(),3);
    vals = {y_[4], (1.0 / L_) , - (R_ / L_) , y_[0]};
    J_.setValues(rcord, ccord, vals);

    
    std::vector<IdxT> ccor2{0, 2, 3, 4};
    std::fill(rcord.begin(), rcord.end(), 4);
    vals = {-y_[3], (1.0 / L_) , -y_[0], - (R_ / L_)};
    J_.setValues(rcord, ccor2, vals);


    //Create -dF/dy'
    std::vector<IdxT> rcordder{3,4};
    std::vector<IdxT> ccordder{3,4};
    std::vector<ScalarT> valsder{-1.0, -1.0};
    COO_Matrix<ScalarT,IdxT> Jacder = COO_Matrix<ScalarT, IdxT>(rcordder, ccordder, valsder,5,5);
    
    //Perform dF/dy + \alpha dF/dy'
    J_.axpy(alpha_, Jacder);

    return 0;
}

template <class ScalarT, typename IdxT>
int MicrogridLoad<ScalarT, IdxT>::evaluateIntegrand()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int MicrogridLoad<ScalarT, IdxT>::initializeAdjoint()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int MicrogridLoad<ScalarT, IdxT>::evaluateAdjointResidual()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int MicrogridLoad<ScalarT, IdxT>::evaluateAdjointIntegrand()
{
    return 0;
}





// Available template instantiations
template class MicrogridLoad<double, long int>;
template class MicrogridLoad<double, size_t>;


} //namespace ModelLib

