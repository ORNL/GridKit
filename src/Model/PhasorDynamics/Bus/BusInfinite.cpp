
#include <iostream>
#include <cmath>

#include <PowerSystemData.hpp>
#include "BusInfinite.hpp"

namespace GridKit
{
namespace PhasorDynamics 
{

/*!
 * @brief Constructor for an infinite (slack) bus.
 *
 * The model is using current balance in Cartesian coordinates.
 * 
 * Arguments to be passed to BusBase:
 * - Number of equations = 0 (size_)
 * - Number of variables = 0 (size_)
 * - Number of quadratures = 0
 * - Number of optimization parameters = 0
 */
template <class ScalarT, typename IdxT>
BusInfinite<ScalarT, IdxT>::BusInfinite()
{
    //std::cout << "Create BusInfinite..." << std::endl;
    //std::cout << "Number of equations is " << size_ << std::endl;

    size_ = 0;
}

/*!
 * @brief BusInfinite constructor.
 *
 * This constructor sets initial values for active and reactive voltage.
 *
 * Arguments to be passed to BusBase:
 * - Number of equations = 0 (size_)
 * - Number of variables = 0 (size_)
 * - Number of quadratures = 0
 * - Number of optimization parameters = 0
 */
template <class ScalarT, typename IdxT>
BusInfinite<ScalarT, IdxT>::BusInfinite(ScalarT Vr, ScalarT Vi)
  : Vr_(Vr), Vi_(Vi)
{
    //std::cout << "Create BusInfinite..." << std::endl;
    //std::cout << "Number of equations is " << size_ << std::endl;

    size_ = 0;
}

/**
 * @brief Construct a new BusInfinite
 * 
 * Arguments to be set in BusBase:
 * - Number of equations = 0 (size_)
 * - Number of variables = 0 (size_)
 * - Number of quadratures = 0
 * - Number of optimization parameters = 0

 * @tparam ScalarT - type of scalar variables
 * @tparam IdxT    - type for vector/matrix indices
 * @param[in] data - structure with bus data 
 */
template <class ScalarT, typename IdxT>
BusInfinite<ScalarT, IdxT>::BusInfinite(BusData& data)
  : BusBase<ScalarT, IdxT>(data.bus_i),
    Vr_(data.Vm * cos(data.Va)),
    Vi_(data.Vm * sin(data.Va))
{
    size_ = 0;
}

template <class ScalarT, typename IdxT>
BusInfinite<ScalarT, IdxT>::~BusInfinite()
{
    //std::cout << "Destroy PQ bus ..." << std::endl;
}

/*!
 * @brief allocate method resizes local solution and residual vectors.
 */
template <class ScalarT, typename IdxT>
int BusInfinite<ScalarT, IdxT>::allocate()
{
    //std::cout << "Nothing to allocate for infinite bus ..." << std::endl;

    return 0;
}


template <class ScalarT, typename IdxT>
int BusInfinite<ScalarT, IdxT>::tagDifferentiable()
{
    return 0;
}


/*!
 * @brief initialize method sets bus variables to stored initial values.
 */
template <class ScalarT, typename IdxT>
int BusInfinite<ScalarT, IdxT>::initialize()
{
    // std::cout << "Initialize BusInfinite..." << std::endl;

    return 0;
}

/*!
 * @brief Reset slack currents to zero.
 *
 * Infinite bus does not compute residuals, so here we just reset
 * current values to zero. Components connected to the infinite bus
 * will add their currents to Ir_ and Ii_. The resultant will be slack
 * current that the infinite bus has to pick up.
 *
 * @warning This implementation assumes bus residuals are always evaluated
 * _before_ component model residuals.
 *
 */
template <class ScalarT, typename IdxT>
int BusInfinite<ScalarT, IdxT>::evaluateResidual()
{
    // std::cout << "Evaluating residual of a PQ bus ...\n";
    Ir_ = 0.0;
    Ii_ = 0.0;
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
int BusInfinite<ScalarT, IdxT>::evaluateJacobian()
{
    return 0;
}

/*!
 * @brief initialize method sets bus variables to stored initial values.
 */
template <class ScalarT, typename IdxT>
int BusInfinite<ScalarT, IdxT>::initializeAdjoint()
{
    // std::cout << "Initialize BusInfinite..." << std::endl;

    return 0;
}

/**
 * @brief BusInfinite only initializes adjoint residual elements to zero.
 * 
 * @tparam ScalarT - data type for the integrand
 * @tparam IdxT    - data type for matrix/vector indices
 * @return int - error code
 */
template <class ScalarT, typename IdxT>
int BusInfinite<ScalarT, IdxT>::evaluateAdjointResidual()
{
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
int BusInfinite<ScalarT, IdxT>::evaluateIntegrand()
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
int BusInfinite<ScalarT, IdxT>::evaluateAdjointIntegrand()
{
    return 0;
}


// Available template instantiations
template class BusInfinite<double, long int>;
template class BusInfinite<double, size_t>;

} // namespace PhasorDynamic
} // namespace GridKit

