
#include <cmath>
#include <iostream>

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusData.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename ScalarT, typename IdxT>
    BusBase<ScalarT, IdxT>::~BusBase() = default;

    template <typename ScalarT, typename IdxT>
    const Model::VariableMonitorBase* BusBase<ScalarT, IdxT>::getMonitor() const
    {
      return monitor_.get();
    }

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
    Bus<ScalarT, IdxT>::Bus(const DataT& data)
      : Vr0_(data.Vr0),
        Vi0_(data.Vi0)
    {
      bus_id_        = data.bus_id;
      size_          = 2;
      monitor_       = std::make_unique<MonitorT>("Bus_" + data.name, data.monitored_variables);
      using Variable = typename DataT::MonitorableVariables;
      monitor_->set(Variable::Vr, [this]
                    { return Vr(); });
      monitor_->set(Variable::Vi, [this]
                    { return Vi(); });
      monitor_->set(Variable::Vm, [this]
                    { return std::sqrt(Vr() * Vr() + Vi() * Vi()); });
      monitor_->set(Variable::Va, [this]
                    { return std::atan2(Vi(), Vr()); });
    }

    template <class ScalarT, typename IdxT>
    Bus<ScalarT, IdxT>::~Bus()
    {
      // std::cout << "Destroy PQ bus ..." << std::endl;
    }

    /*!
     * @brief allocate method resizes local solution and residual vectors.
     */
    template <class ScalarT, typename IdxT>
    int Bus<ScalarT, IdxT>::allocate()
    {
      // Temporary while we use std::vector in the code
      size_t size = static_cast<size_t>(size_);

      // Resize component model data
      f_.resize(size);
      y_.resize(size);
      yp_.resize(size);
      tag_.resize(size);
      variable_indices_.resize(size);
      residual_indices_.resize(size);

      // Default variable and residual index mapping to local index
      for (IdxT j = 0; j < size_; ++j)
      {
        this->setVariableIndex(j, j);
        this->setResidualIndex(j, j);
      }

      return 0;
    }

    /**
     * @brief Set the bus ID
     */
    template <class ScalarT, typename IdxT>
    int Bus<ScalarT, IdxT>::setBusID(IdxT bus_id)
    {
      bus_id_ = bus_id;
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
      y_[0]  = Vr0_;
      y_[1]  = Vi0_;
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
  } // namespace PhasorDynamics
} // namespace GridKit
