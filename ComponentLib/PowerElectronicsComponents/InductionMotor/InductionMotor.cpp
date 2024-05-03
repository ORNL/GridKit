


#include <iostream>
#include <cmath>
#include <vector>
#include "InductionMotor.hpp"


namespace ModelLib {



/*!
 * @brief Constructor for a constant InductionMotor model
 *
 * Calls default ModelEvaluatorImpl constructor.
 * @todo create a test case utilizing the component.
 * @todo create a unit test to check correctness of component
 */

template <class ScalarT, typename IdxT>
InductionMotor<ScalarT, IdxT>::InductionMotor(IdxT id, ScalarT Lls, ScalarT Rs, ScalarT Llr, ScalarT Rr, ScalarT Lms, ScalarT RJ, ScalarT P)
  : Lls_(Lls),
    Rs_(Rs),
    Llr_(Llr),
    Rr_(Rr),
    Lms_(Lms),
    RJ_(RJ),
    P_(P)
{
    size_ = 10;
    n_intern_ = 5;
    n_extern_ = 5;
    extern_indices_ = {0,1,2,3,4};
    idc_ = id;
}

template <class ScalarT, typename IdxT>
InductionMotor<ScalarT, IdxT>::~InductionMotor()
{
}

/*!
 * @brief allocate method computes sparsity pattern of the Jacobian.
 */
template <class ScalarT, typename IdxT>
int InductionMotor<ScalarT, IdxT>::allocate()
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
int InductionMotor<ScalarT, IdxT>::initialize()
{
    return 0;
}

/*
 * \brief Identify differential variables
 */
template <class ScalarT, typename IdxT>
int InductionMotor<ScalarT, IdxT>::tagDifferentiable()
{
    return 0;
}

/**
 * @brief Contributes to the resisdual
 *
 */
template <class ScalarT, typename IdxT>
int InductionMotor<ScalarT, IdxT>::evaluateResidual()
{
    
    f_[0] = y_[5] + y_[7];
    f_[1] = (-1.0/2.0) * y_[5] - (sqrt(3.0)/2.0)*y_[6] + y_[7];
    f_[2] = (-1.0/2.0) * y_[5] + (sqrt(3.0)/2.0)*y_[6] + y_[7];
    f_[3] = RJ_ * yp_[3] - (3.0/4.0)*P_*Lms_ * (y_[5]*y_[9] - y_[6]*y_[8]);
    f_[4] = yp_[4] - y_[3];
    f_[5] = (1.0/3.0)*(2.0* y_[0] - y_[1] - y_[2]) - Rs_*y_[5] - (Lls_ + Lms_) * yp_[5] - Lms_ * yp_[6];
    f_[6] = (1.0/sqrt(3.0))*(-y_[1] + y_[2]) - Rs_*y_[6] - (Lls_ + Lms_) * yp_[6] - Lms_ * yp_[5];
    f_[7] = (y_[0] + y_[1] + y_[2])/3.0 - Rs_*y_[7] - Lls_ * yp_[7];
    f_[8] = Rr_*y_[8] + (Llr_ + Lms_)*yp_[8] + Lms_ * yp_[5] - (P_/2)*y_[3]*((Llr_+Lms_)*y_[9] + Lms_*y_[6]);
    f_[9] = Rr_*y_[9] + (Llr_ + Lms_)*yp_[9] + Lms_ * yp_[6] + (P_/2)*y_[3]*((Llr_+Lms_)*y_[8] + Lms_*y_[5]);
    return 0;
}

/**
 * @brief Compute component Jacobian
 * 
 * @todo need to implement
 * 
 * @tparam ScalarT 
 * @tparam IdxT 
 * @return int 
 */
template <class ScalarT, typename IdxT>
int InductionMotor<ScalarT, IdxT>::evaluateJacobian()
{
    
    return 0;
}

template <class ScalarT, typename IdxT>
int InductionMotor<ScalarT, IdxT>::evaluateIntegrand()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int InductionMotor<ScalarT, IdxT>::initializeAdjoint()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int InductionMotor<ScalarT, IdxT>::evaluateAdjointResidual()
{
    return 0;
}

template <class ScalarT, typename IdxT>
int InductionMotor<ScalarT, IdxT>::evaluateAdjointIntegrand()
{
    return 0;
}





// Available template instantiations
template class InductionMotor<double, long int>;
template class InductionMotor<double, size_t>;


} //namespace ModelLib

