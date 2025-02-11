
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
 * @brief PQ bus does not compute residuals, so here we just reset residual values.
 *
 * @warning This implementation assumes bus residuals are always evaluated
 * _before_ component model residuals.
 *
 */
template <class ScalarT, typename IdxT>
int BusInfinite<ScalarT, IdxT>::evaluateResidual()
{
    // std::cout << "Evaluating residual of a PQ bus ...\n";
    f_[0] = 0.0;
    f_[1] = 0.0;
    return 0;
}


/*!
 * @brief initialize method sets bus variables to stored initial values.
 */
template <class ScalarT, typename IdxT>
int BusInfinite<ScalarT, IdxT>::initializeAdjoint()
{
    // std::cout << "Initialize BusInfinite..." << std::endl;
    yB_[0] = 0.0;
    yB_[1] = 0.0;
    ypB_[0] = 0.0;
    ypB_[1] = 0.0;

    return 0;
}

template <class ScalarT, typename IdxT>
int BusInfinite<ScalarT, IdxT>::evaluateAdjointResidual()
{
    fB_[0] = 0.0;
    fB_[1] = 0.0;

    return 0;
}

// Available template instantiations
template class BusInfinite<double, long int>;
template class BusInfinite<double, size_t>;

} // namespace PhasorDynamic
} // namespace GridKit

