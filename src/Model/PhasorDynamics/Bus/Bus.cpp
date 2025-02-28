
#include <iostream>
#include <cmath>

#include <PowerSystemData.hpp>
#include "Bus.hpp"

namespace GridKit
{
namespace PhasorDynamics 
{

/*!
 * @brief Constructor for a phasor dynamics bus.
 *
 * The model is using current balance in Cartesian coordinates.
 * 
 * @todo Arguments that should be passed to ModelEvaluatorImpl constructor:
 * - Number of equations = 2 (size_)
 * - Number of variables = 2 (size_)
 * - Number of quadratures = 0
 * - Number of optimization parameters = 0
 */
template <class ScalarT, typename IdxT>
Bus<ScalarT, IdxT>::Bus()
  : Vr0_(0.0), Vi0_(0.0)
{
    //std::cout << "Create Bus..." << std::endl;
    //std::cout << "Number of equations is " << size_ << std::endl;

    size_ = 2;
}

/*!
 * @brief Bus constructor.
 *
 * This constructor sets initial values for active and reactive voltage.
 *
 * @todo Arguments that should be passed to ModelEvaluatorImpl constructor:
 * - Number of equations = 2 (size_)
 * - Number of variables = 2 (size_)
 * - Number of quadratures = 0
 * - Number of optimization parameters = 0
 */
template <class ScalarT, typename IdxT>
Bus<ScalarT, IdxT>::Bus(ScalarT Vr, ScalarT Vi)
  : Vr0_(Vr), Vi0_(Vi)
{
    //std::cout << "Create Bus..." << std::endl;
    //std::cout << "Number of equations is " << size_ << std::endl;

    size_ = 2;
}

/**
 * @brief Construct a new Bus
 * 
 * @tparam ScalarT - type of scalar variables
 * @tparam IdxT    - type for vector/matrix indices
 * @param[in] data - structure with bus data 
 */
template <class ScalarT, typename IdxT>
Bus<ScalarT, IdxT>::Bus(BusData& data)
  : BusBase<ScalarT, IdxT>(data.bus_i),
    Vr0_(data.Vm * cos(data.Va)),
    Vi0_(data.Vm * sin(data.Va))
{
    //std::cout << "Create Bus..." << std::endl;
    //std::cout << "Number of equations is " << size_ << std::endl;

    size_ = 2;
}

template <class ScalarT, typename IdxT>
Bus<ScalarT, IdxT>::~Bus()
{
    //std::cout << "Destroy PQ bus ..." << std::endl;
}

/*!
 * @brief allocate method resizes local solution and residual vectors.
 */
template <class ScalarT, typename IdxT>
int Bus<ScalarT, IdxT>::allocate()
{
    //std::cout << "Allocate PQ bus ..." << std::endl;
    f_.resize(size_);
    y_.resize(size_);
    yp_.resize(size_);
    tag_.resize(size_);

    fB_.resize(size_);
    yB_.resize(size_);
    ypB_.resize(size_);

    return 0;
}


/*!
 * @brief Bus variables are algebraic.
 */
template <class ScalarT, typename IdxT>
int Bus<ScalarT, IdxT>::tagDifferentiable()
{
    tag_[0] = false;
    tag_[1] = false;
    return 0;
}


/*!
 * @brief initialize method sets bus variables to stored initial values.
 */
template <class ScalarT, typename IdxT>
int Bus<ScalarT, IdxT>::initialize()
{
    // std::cout << "Initialize Bus..." << std::endl;
    y_[0] = Vr0_;
    y_[1] = Vi0_;
    yp_[0] = 0.0;
    yp_[1] = 0.0;

    return 0;
}

/*!
 * @brief PQ bus does not compute residuals, so here we just reset residual values.
 *
 * @warning This implementation assumes bus residuals are always evaluated
 * _before_ component model residuals.
 *
 */
template <class ScalarT, typename IdxT>
int Bus<ScalarT, IdxT>::evaluateResidual()
{
    // std::cout << "Evaluating residual of a PQ bus ...\n";
    f_[0] = 0.0;
    f_[1] = 0.0;
    return 0;
}

/**
 * @brief Jacobian evaluation not implemented
 * 
 * @tparam ScalarT - data type for Jacobian elements
 * @tparam IdxT    - data type for matrix indices
 * @return int - error code
 */
template <class ScalarT, typename IdxT>
int Bus<ScalarT, IdxT>::evaluateJacobian()
{
    return 0;
}


/*!
 * @brief initialize method sets bus variables to stored initial values.
 */
template <class ScalarT, typename IdxT>
int Bus<ScalarT, IdxT>::initializeAdjoint()
{
    // std::cout << "Initialize Bus..." << std::endl;
    yB_[0] = 0.0;
    yB_[1] = 0.0;
    ypB_[0] = 0.0;
    ypB_[1] = 0.0;

    return 0;
}

/**
 * @brief Bus only initializes adjoint residual elements to zero.
 * 
 * @tparam ScalarT - data type for the integrand
 * @tparam IdxT    - data type for matrix/vector indices
 * @return int - error code
 */
template <class ScalarT, typename IdxT>
int Bus<ScalarT, IdxT>::evaluateAdjointResidual()
{
    fB_[0] = 0.0;
    fB_[1] = 0.0;

    return 0;
}

/**
 * @brief Quadrature evaluation not implemented
 * 
 * @tparam ScalarT - data type for the integrand
 * @tparam IdxT    - data type for matrix/vector indices
 * @return int - error code
 */
template <class ScalarT, typename IdxT>
int Bus<ScalarT, IdxT>::evaluateIntegrand()
{
    return 0;
}

/**
 * @brief Adjoint quadrature evaluation not implemented
 * 
 * @tparam ScalarT - data type for the integrand
 * @tparam IdxT    - data type for matrix/vector indices
 * @return int - error code
 */
template <class ScalarT, typename IdxT>
int Bus<ScalarT, IdxT>::evaluateAdjointIntegrand()
{
    return 0;
}

// Available template instantiations
template class Bus<double, long int>;
template class Bus<double, size_t>;

} // namespace PhasorDynamic
} // namespace GridKit

