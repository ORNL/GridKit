
#include <cmath>
#include <iostream>

#include <GridKit/Model/PhasorDynamics/Bus/BusData.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <GridKit/Model/PhasorDynamics/BusBaseImpl.hpp>
#include <GridKit/Model/PhasorDynamics/GridElementImpl.hpp>

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
    template <typename ScalarP, typename IdxP>
    BusInfinite<ScalarP, IdxP>::BusInfinite()
    {
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
    template <typename ScalarP, typename IdxP>
    BusInfinite<ScalarP, IdxP>::BusInfinite(ScalarT Vr, ScalarT Vi)
      : Vr_(Vr), Vi_(Vi)
    {
    }

    /**
     * @brief Construct a new BusInfinite
     *
     * Arguments to be set in BusBase:
     * - Number of equations = 0 (size_)
     * - Number of variables = 0 (size_)

     * @param[in] data - structure with bus data
     */
    template <typename ScalarP, typename IdxP>
    BusInfinite<ScalarP, IdxP>::BusInfinite(const ModelDataT& data)
      : ConnectedElement<BusInfinite>(data),
        Vr_(data.Vr0),
        Vi_(data.Vi0)
    {
      using Variable = typename BusData<RealT, IdxT>::MonitorableVariables;

      monitor_->set(Variable::Vr, [this]
                    { return Vr(); });
      monitor_->set(Variable::Vi, [this]
                    { return Vi(); });
      monitor_->set(Variable::Vm, [this]
                    { return std::sqrt(Vr() * Vr() + Vi() * Vi()); });
      monitor_->set(Variable::Va, [this]
                    { return std::atan2(Vi(), Vr()); });
    }

    template <typename ScalarP, typename IdxP>
    BusInfinite<ScalarP, IdxP>::~BusInfinite()
    {
    }

    /**
     * @brief Set the bus ID
     */
    template <typename ScalarP, typename IdxP>
    int BusInfinite<ScalarP, IdxP>::setBusID(IdxT bus_id)
    {
      bus_id_ = bus_id;
      return 0;
    }

    /*!
     * @brief allocate method resizes local solution and residual vectors.
     */
    template <typename ScalarP, typename IdxP>
    int BusInfinite<ScalarP, IdxP>::allocate()
    {
      return 0;
    }

    /**
     * @brief Tag differentiable variables
     */
    template <typename ScalarP, typename IdxP>
    int BusInfinite<ScalarP, IdxP>::tagDifferentiable()
    {
      return 0;
    }

    /*!
     * @brief initialize method sets bus variables to stored initial values.
     */
    template <typename ScalarP, typename IdxP>
    int BusInfinite<ScalarP, IdxP>::initialize()
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
    template <typename ScalarP, typename IdxP>
    int BusInfinite<ScalarP, IdxP>::evaluateResidual()
    {
      Ir_ = 0.0;
      Ii_ = 0.0;
      return 0;
    }

    /**
     * @brief There is no Jacobian for slack variables
     *
     * @return int - error code
     */
    template <typename ScalarP, typename IdxP>
    int BusInfinite<ScalarP, IdxP>::evaluateJacobian()
    {
      return 0;
    }

  } // namespace PhasorDynamics
} // namespace GridKit
