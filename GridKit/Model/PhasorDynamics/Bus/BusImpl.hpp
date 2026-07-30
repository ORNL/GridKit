
#include <cmath>
#include <iostream>

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusData.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename scalar_type, typename index_type>
    BusBase<scalar_type, index_type>::~BusBase() = default;

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* BusBase<scalar_type, index_type>::getMonitor() const
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
    template <typename scalar_type, typename index_type>
    Bus<scalar_type, index_type>::Bus()
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
    template <typename scalar_type, typename index_type>
    Bus<scalar_type, index_type>::Bus(ScalarT Vr, ScalarT Vi)
      : Vr0_(Vr), Vi0_(Vi)
    {
      size_ = 2;
    }

    /**
     * @brief Construct a new Bus
     *
     * @param[in] data - structure with bus data
     */
    template <typename scalar_type, typename index_type>
    Bus<scalar_type, index_type>::Bus(const ModelDataT& data)
      : Vr0_(data.Vr0),
        Vi0_(data.Vi0)
    {
      bus_id_        = data.bus_id;
      size_          = 2;
      monitor_       = std::make_unique<MonitorT>("Bus_" + data.name, data.monitored_variables);
      using Variable = typename ModelDataT::MonitorableVariables;
      monitor_->set(Variable::Vr, [this]
                    { return Vr(); });
      monitor_->set(Variable::Vi, [this]
                    { return Vi(); });
      monitor_->set(Variable::Vm, [this]
                    { return std::sqrt(Vr() * Vr() + Vi() * Vi()); });
      monitor_->set(Variable::Va, [this]
                    { return std::atan2(Vi(), Vr()); });
    }

    template <typename scalar_type, typename index_type>
    Bus<scalar_type, index_type>::~Bus() = default;

    /*!
     * @brief Allocate bus storage and index maps.
     */
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::allocate()
    {
      if (!allocated_)
      {
        this->allocateVectors(size_);
      }
      size_t size = static_cast<size_t>(size_);

      tag_.resize(size);

      variable_indices_.resize(size);
      residual_indices_.resize(size);

      // Default variable and residual index mapping to local index
      for (IdxT j = 0; j < size_; ++j)
      {
        this->setVariableIndex(j, j);
        this->setResidualIndex(j, j);
      }

      // Publish voltage signals assigned to this bus
      static constexpr auto VR = BusInternalVariables::VR;
      static constexpr auto VI = BusInternalVariables::VI;

      if (signals_.template isAssigned<VR>())
      {
        signals_.template getSignalNode<VR>()->set(
            &y_.getData()[0], &this->getVariableIndex(0), &this->getResidualIndex(0));
      }
      if (signals_.template isAssigned<VI>())
      {
        signals_.template getSignalNode<VI>()->set(
            &y_.getData()[1], &this->getVariableIndex(1), &this->getResidualIndex(1));
      }

      allocated_ = true;
      return 0;
    }

    /**
     * @brief Set the bus ID
     */
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::setBusID(IdxT bus_id)
    {
      bus_id_ = bus_id;
      return 0;
    }

    /*!
     * @brief Bus variables are algebraic.
     */
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::tagDifferentiable()
    {
      tag_[0] = false;
      tag_[1] = false;
      return 0;
    }

    /**
     * @brief Compute the absolute tolerance for each variable in the model
     *
     * @param rel_tol The relative tolerance which can be used to pick the
     *        absolute tolerance.
     * @tparam scalar_type Scalar data type
     * @tparam index_type Index data type
     * @return int 0 if successful, non-zero otherwise.
     *
     * This represents a "noise" level close to zero for which pure relative
     * error cannot be used.
     */
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return 0;
    }

    /*!
     * @brief initialize method sets bus variables to stored initial values.
     */
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::initialize()
    {
      // std::cout << "Initialize Bus..." << std::endl;
      auto* y  = y_.getData();
      auto* yp = yp_.getData();

      y[0]  = Vr0_;
      y[1]  = Vi0_;
      yp[0] = 0.0;
      yp[1] = 0.0;

      y_.setDataUpdated();
      yp_.setDataUpdated();

      return 0;
    }

    /*!
     * @brief Assign the bus current balance rows.
     *
     * The bus owns the KCL rows and assigns them to zero; connected
     * components accumulate their current contributions in the external
     * residual phase.
     */
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::evaluateInternalResidual()
    {
      auto* f = f_.getData();

      f[0] = 0.0;
      f[1] = 0.0;
      f_.setDataUpdated();
      return 0;
    }

    /*!
     * @brief Evaluate the internal residual and external residual
     * contributions.
     */
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::evaluateResidual()
    {
      evaluateInternalResidual();
      this->evaluateExternalResidual();

      return 0;
    }
  } // namespace PhasorDynamics
} // namespace GridKit
