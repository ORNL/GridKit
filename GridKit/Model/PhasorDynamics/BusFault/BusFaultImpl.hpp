#pragma once

#include <GridKit/Constants.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/BusFault/BusFault.hpp>
#include <GridKit/Model/PhasorDynamics/BusFault/BusFaultData.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Constructor for a bus fault
     *
     * Model size:
     * - Number of equations = 0
     * - Number of internal variables = 0
     */
    template <typename scalar_type, typename index_type>
    BusFault<scalar_type, index_type>::BusFault(BusT* bus)
      : bus_(bus), R_(0), X_(0.01), bus_id_(0)
    {
      size_ = 0;
      setDerivedParams();
    }

    /**
     * @brief Construct a new BusFault from model data
     *
     * @param bus - pointer to the faulted bus
     * @param data - bus fault model data
     */
    template <typename scalar_type, typename index_type>
    BusFault<scalar_type, index_type>::BusFault(BusT* bus, const ModelDataT& data)
      : bus_(bus),
        monitor_(std::make_unique<MonitorT>(data))
    {
      using Parameter = typename ModelDataT::Parameters;
      using Buses     = typename ModelDataT::Buses;

      if (data.parameters.contains(Parameter::R))
      {
        R_ = std::get<RealT>(data.parameters.at(Parameter::R));
      }

      if (data.parameters.contains(Parameter::X))
      {
        X_ = std::get<RealT>(data.parameters.at(Parameter::X));
      }

      if (data.buses.contains(Buses::bus))
      {
        bus_id_ = data.buses.at(Buses::bus);
      }

      initializeMonitor();

      size_ = 0;
      setDerivedParams();
    }

    template <typename scalar_type, typename index_type>
    BusFault<scalar_type, index_type>::~BusFault()
    {
    }

    /**
     * @brief Set the component ID
     */
    template <typename scalar_type, typename index_type>
    int BusFault<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /*!
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <typename scalar_type, typename index_type>
    int BusFault<scalar_type, index_type>::allocate()
    {
      if (!allocated_)
      {
        this->allocateVectors(size_);
      }
      auto size = static_cast<std::size_t>(size_);

      variable_indices_.resize(size);
      residual_indices_.resize(size);

      // Resize coupling data
      wb_.resize(2);
      h_.resize(2);

      allocated_ = true;
      return 0;
    }

    /**
     * Initialization of the bus fault model
     *
     */
    template <typename scalar_type, typename index_type>
    int BusFault<scalar_type, index_type>::initialize()
    {
      status_ = ZERO<RealT>;
      if (signals_.template isAttached<BusFaultExternalVariables::STATUS>())
      {
        status_ = signals_.template readExternalVariable<BusFaultExternalVariables::STATUS>();
      }

      return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <typename scalar_type, typename index_type>
    int BusFault<scalar_type, index_type>::tagDifferentiable()
    {
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
    int BusFault<scalar_type, index_type>::setAbsoluteTolerance(RealT)
    {
      return 0;
    }

    /**
     * @brief Bus residual contribution from bus variables
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline int BusFault<scalar_type, index_type>::evaluateBusResidual(
        [[maybe_unused]] const ScalarT* y,
        [[maybe_unused]] const ScalarT* yp,
        const ScalarT*                  wb,
        ScalarT*                        h)
    {
      const ScalarT Vr = wb[0];
      const ScalarT Vi = wb[1];

      h[0] = g_ * Vr - b_ * Vi;
      h[1] = b_ * Vr + g_ * Vi;

      return 0;
    }

    /**
     * @brief Residual contribution of the fault is computed and pushed to the faulted bus.
     *
     */
    template <typename scalar_type, typename index_type>
    int BusFault<scalar_type, index_type>::evaluateResidual()
    {
      status_ = ZERO<RealT>;
      if (signals_.template isAttached<BusFaultExternalVariables::STATUS>())
      {
        status_ = signals_.template readExternalVariable<BusFaultExternalVariables::STATUS>();
      }

      ScalarT ir{0.0};
      ScalarT ii{0.0};

      faultCurrent(ir, ii);

      Ir() += ir;
      Ii() += ii;

      if (bus_->size() > 0)
      {
        bus_->getResidual().setDataUpdated();
      }

      return 0;
    }

    template <typename scalar_type, typename index_type>
    void BusFault<scalar_type, index_type>::faultCurrent(ScalarT& Ir, ScalarT& Ii)
    {
      const ScalarT wb[2]{Vr(), Vi()};
      ScalarT       h[2];
      evaluateBusResidual(nullptr, nullptr, wb, h);

      Ir = status_ * h[0];
      Ii = status_ * h[1];
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* BusFault<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    template <typename scalar_type, typename index_type>
    void BusFault<scalar_type, index_type>::initializeMonitor()
    {
      using Variable = typename ModelDataT::MonitorableVariables;

      monitor_->set(Variable::status, [this]
                    { return status_; });
      monitor_->set(Variable::ir, [this]
                    {
                      ScalarT ir{0.0};
                      ScalarT ii{0.0};
                      faultCurrent(ir, ii);
                      return ir; });
      monitor_->set(Variable::ii, [this]
                    {
                      ScalarT ir{0.0};
                      ScalarT ii{0.0};
                      faultCurrent(ir, ii);
                      return ii; });
    }

    /**
     * @brief Derived parameters
     *
     */
    template <typename scalar_type, typename index_type>
    void BusFault<scalar_type, index_type>::setDerivedParams()
    {
      g_ = ZERO<RealT>;
      b_ = ZERO<RealT>;

      const RealT denom = R_ * R_ + X_ * X_;
      if (denom == ZERO<RealT>)
      {
        return;
      }

      g_ = -R_ / denom;
      b_ = X_ / denom;
    }

  } // namespace PhasorDynamics
} // namespace GridKit
