
#include <iostream>
#include <cmath>
#include <vector>
#include "MicrogridBusDQ.hpp"

namespace ModelLib {

/*!
 * @brief Constructor for a constant MicrogridBusDQ model
 *
 * Calls default ModelEvaluatorImpl constructor.
 * 
 * In DQ space
 * Each microgrid line has a virtual resistance RN
 * Model is from paper: "Modeling, Analysis and Testing of Autonomous Operation
 * of an Inverter-Based Microgrid", Nagaraju Pogaku, Milan Prodanovic, and
 * Timothy C. Green, Section E
 */
template <class ScalarT, typename IdxT>
MicrogridBusDQ<ScalarT, IdxT>::MicrogridBusDQ(IdxT id, ScalarT RN)
  : RN_(RN) 
{
    // externals [vbus_d, vbus_q]
    size_ = 2;
    n_intern_ = 0;
    n_extern_ = 2;
    extern_indices_ = {0,1};
    idc_ = id;
}

template <class ScalarT, typename IdxT>
MicrogridBusDQ<ScalarT, IdxT>::~MicrogridBusDQ()
{
}

/*!
 * @brief allocate method computes sparsity pattern of the Jacobian.
 */
template <class ScalarT, typename IdxT>
int MicrogridBusDQ<ScalarT, IdxT>::allocate()
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
int MicrogridBusDQ<ScalarT, IdxT>::initialize()
{
    return 0;
}

/*
 * \brief Identify differential variables
 */
template <class ScalarT, typename IdxT>
int MicrogridBusDQ<ScalarT, IdxT>::tagDifferentiable()
{
    return 0;
}

/**
 * @brief Evaluate residual of microgrid line
 * This model has "Virtual resistors". The voltage of the bus divided by its virtual resistance.
 * The components are external to allow for outside components to add inductances to the terms.
 *
 * refernce to equations in class header
 * 
 */
template <class ScalarT, typename IdxT>
int MicrogridBusDQ<ScalarT, IdxT>::evaluateResidual()
{
    //bus voltage
    f_[0] = -y_[0] / RN_;
    f_[1] = -y_[1] / RN_;

    return 0;
}

/**
 * @brief Generate Jacobian
 * 
 * @tparam ScalarT 
 * @tparam IdxT 
 * @return int 
 */
template <class ScalarT, typename IdxT>
int MicrogridBusDQ<ScalarT, IdxT>::evaluateJacobian()
{
    J_.zeroMatrix();

    //Create dF/dy
    std::vector<IdxT> rtemp{0,1};
    std::vector<IdxT> ctemp{0,1};
    std::vector<ScalarT> vals{-1.0 / RN_,-1.0 / RN_};
    J_.setValues(rtemp, ctemp, vals);

    return 0;
}

template <class ScalarT, typename IdxT>
int MicrogridBusDQ<ScalarT, IdxT>::evaluateIntegrand()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int MicrogridBusDQ<ScalarT, IdxT>::initializeAdjoint()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int MicrogridBusDQ<ScalarT, IdxT>::evaluateAdjointResidual()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int MicrogridBusDQ<ScalarT, IdxT>::evaluateAdjointIntegrand()
{
    return 0;
}


// Available template instantiations
template class MicrogridBusDQ<double, long int>;
template class MicrogridBusDQ<double, size_t>;


} //namespace ModelLib

