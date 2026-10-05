#pragma once

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Load/LoadZIP/LoadZIP.hpp>
#include <GridKit/Model/PhasorDynamics/Load/LoadZIP/LoadZIPData.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Constructor for a zip load
     *
     * Model size:
     * - Number of equations = 0
     * - Number of internal variables = 0
     */
    template <typename scalar_type, typename index_type>
    LoadZIP<scalar_type, index_type>::LoadZIP(BusT* bus)
      : bus_(bus)
    {
      size_ = 0;
      setDerivedParams();
    }

    template <typename scalar_type, typename index_type>
    LoadZIP<scalar_type, index_type>::LoadZIP(BusT* bus, RealT Pnom, RealT Qnom, RealT alphaI, RealT alphaP)
      : bus_(bus),
        Pnom_(Pnom),
        Qnom_(Qnom),
        alphaI_(alphaI),
        alphaP_(alphaP)
    {
      size_ = 0;
      setDerivedParams();
    }

    template <typename scalar_type, typename index_type>
    LoadZIP<scalar_type, index_type>::LoadZIP(BusT*             bus,
                                              const ModelDataT& data)
      : bus_(bus),
        monitor_(std::make_unique<MonitorT>(data))
    {
      initializeParameters(data);
      initializeMonitor();
      size_ = 0;
    }

    template <typename scalar_type, typename index_type>
    LoadZIP<scalar_type, index_type>::~LoadZIP()
    {
    }

    template <typename scalar_type, typename index_type>
    void LoadZIP<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Parameter = typename ModelDataT::Parameters;
      if (data.parameters.contains(Parameter::Pnom))
      {
        Pnom_ = std::get<RealT>(data.parameters.at(Parameter::Pnom));
      }

      if (data.parameters.contains(Parameter::Qnom))
      {
        Qnom_ = std::get<RealT>(data.parameters.at(Parameter::Qnom));
      }

      if (data.parameters.contains(Parameter::alphaI))
      {
        alphaI_ = std::get<RealT>(data.parameters.at(Parameter::alphaI));
      }

      if (data.parameters.contains(Parameter::alphaP))
      {
        alphaP_ = std::get<RealT>(data.parameters.at(Parameter::alphaP));
      }

      setDerivedParams();
    }

    /**
     * @brief Set the component ID
     */
    template <typename scalar_type, typename index_type>
    int LoadZIP<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /*!
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <typename scalar_type, typename index_type>
    int LoadZIP<scalar_type, index_type>::allocate()
    {
      if (!allocated_)
      {
        this->allocateVectors(size_);
      }

      auto size = static_cast<size_t>(size_); // avoid compiler warnings

      tag_.resize(size);

      variable_indices_.resize(size);
      residual_indices_.resize(size);

      // Default variable and residual index mapping to local index
      for (IdxT j = 0; j < size_; ++j)
      {
        this->setVariableIndex(j, j);
        this->setResidualIndex(j, j);
      }

      // Resize coupling data
      wb_.resize(2);
      h_.resize(2);

      allocated_ = true;
      return 0;
    }

    /**
     * Initialization of the load model
     *
     */
    template <typename scalar_type, typename index_type>
    int LoadZIP<scalar_type, index_type>::initialize()
    {
      const ScalarT vr = Vr();
      const ScalarT vi = Vi();

      // The nominal dispatch, initialized bus voltage, and ZIP anchor are one
      // initial condition, so derive them together on every reset. Anchoring
      // at the initialized voltage makes the dispatch exact there for any
      // load fractions.
      const RealT vm0 = static_cast<RealT>(std::sqrt(vr * vr + vi * vi));
      if (!(vm0 > RealT{0}) || !std::isfinite(vm0))
      {
        return 1;
      }
      Vnom_ = vm0;
      setDerivedParams();

      return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <typename scalar_type, typename index_type>
    int LoadZIP<scalar_type, index_type>::tagDifferentiable()
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
    int LoadZIP<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return 0;
    }

    /**
     * @brief Bus residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int LoadZIP<scalar_type, index_type>::evaluateBusResidual(
        [[maybe_unused]] const ScalarT* y,
        [[maybe_unused]] const ScalarT* yp,
        const ScalarT*                  wb,
        ScalarT*                        h)
    {
      const ScalarT Vr    = wb[0];
      const ScalarT Vi    = wb[1];
      const RealT   Vnom2 = Vnom_ * Vnom_;
      const ScalarT V2    = Vr * Vr + Vi * Vi;
      const ScalarT V     = std::sqrt(V2);
      const ScalarT zip   = alphaZ_ + alphaI_ * Vnom_ / V + alphaP_ * Vnom2 / V2;

      h[0] = -(G_ * Vr + B_ * Vi) * zip;
      h[1] = -(G_ * Vi - B_ * Vr) * zip;

      return 0;
    }

    /**
     * @brief Residual contribution of the load is pushed to the bus.
     *
     */
    template <typename scalar_type, typename index_type>
    int LoadZIP<scalar_type, index_type>::evaluateResidual()
    {
      auto* wb = wb_.getData();
      wb[0]    = Vr();
      wb[1]    = Vi();

      auto* h = h_.getData();
      evaluateBusResidual(y_.getData(), yp_.getData(), wb, h);
      Ir() += h[0];
      Ii() += h[1];
      if (bus_->size() > 0)
      {
        bus_->getResidual().setDataUpdated();
      }

      return 0;
    }

    /**
     * @brief Terminal current injected at the present bus voltage
     *
     */
    template <typename scalar_type, typename index_type>
    void LoadZIP<scalar_type, index_type>::terminalCurrent(ScalarT& ir, ScalarT& ii)
    {
      const ScalarT wb[2] = {Vr(), Vi()};
      ScalarT       h[2];
      evaluateBusResidual(y_.getData(), yp_.getData(), wb, h);
      ir = h[0];
      ii = h[1];
    }

    /**
     * @brief Derived parameters
     *
     */
    template <typename scalar_type, typename index_type>
    void LoadZIP<scalar_type, index_type>::setDerivedParams()
    {
      const RealT Vnom2 = Vnom_ * Vnom_;

      G_      = Pnom_ / Vnom2;
      B_      = Qnom_ / Vnom2;
      alphaZ_ = ONE<RealT> - alphaI_ - alphaP_;
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* LoadZIP<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    template <typename scalar_type, typename index_type>
    void LoadZIP<scalar_type, index_type>::initializeMonitor()
    {
      using Variable = typename ModelDataT::MonitorableVariables;

      monitor_->set(Variable::ir, [this]
                    {
                      ScalarT ir;
                      ScalarT ii;
                      terminalCurrent(ir, ii);
                      return ir; });
      monitor_->set(Variable::ii, [this]
                    {
                      ScalarT ir;
                      ScalarT ii;
                      terminalCurrent(ir, ii);
                      return ii; });
      monitor_->set(Variable::im, [this]
                    {
                      ScalarT ir;
                      ScalarT ii;
                      terminalCurrent(ir, ii);
                      return std::sqrt(ir * ir + ii * ii); });
      monitor_->set(Variable::p, [this]
                    {
                      ScalarT ir;
                      ScalarT ii;
                      terminalCurrent(ir, ii);
                      return Vr() * ir + Vi() * ii; });
      monitor_->set(Variable::q, [this]
                    {
                      ScalarT ir;
                      ScalarT ii;
                      terminalCurrent(ir, ii);
                      return Vi() * ir - Vr() * ii; });
    }

  } // namespace PhasorDynamics
} // namespace GridKit
