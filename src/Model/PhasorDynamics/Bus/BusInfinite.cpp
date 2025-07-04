
#include "BusInfinite.hpp"

#include <cmath>
#include <iostream>

#include <Model/PhasorDynamics/Bus/BusData.hpp>

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
     */
    template <class ScalarT, typename IdxT>
    BusInfinite<ScalarT, IdxT>::BusInfinite()
    {
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
     */
    template <class ScalarT, typename IdxT>
    BusInfinite<ScalarT, IdxT>::BusInfinite(ScalarT Vr, ScalarT Vi)
      : Vr_(Vr), Vi_(Vi)
    {
      size_ = 0;
    }

    /**
     * @brief Construct a new BusInfinite
     *
     * Arguments to be set in BusBase:
     * - Number of equations = 0 (size_)
     * - Number of variables = 0 (size_)

     * @tparam ScalarT - type of scalar variables
     * @tparam IdxT    - type for vector/matrix indices
     * @param[in] data - structure with bus data
     */
    template <class ScalarT, typename IdxT>
    BusInfinite<ScalarT, IdxT>::BusInfinite(const DataT& data)
      : BusBase<ScalarT, IdxT>(data.bus_id),
        Vr_(data.Vr0),
        Vi_(data.Vi0)
    {
      size_ = 0;
    }

    template <class ScalarT, typename IdxT>
    BusInfinite<ScalarT, IdxT>::~BusInfinite()
    {
    }

    /*!
     * @brief allocate method resizes local solution and residual vectors.
     */
    template <class ScalarT, typename IdxT>
    int BusInfinite<ScalarT, IdxT>::allocate()
    {
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

    // Available template instantiations
    template class BusInfinite<double, long int>;
    template class BusInfinite<double, size_t>;
    template class BusInfinite<DependencyTracking::Variable, long int>;
    template class BusInfinite<DependencyTracking::Variable, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
