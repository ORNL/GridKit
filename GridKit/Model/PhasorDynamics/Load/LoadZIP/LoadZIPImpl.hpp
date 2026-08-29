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

      variable_indices_.resize(size);
      residual_indices_.resize(size);

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
    int LoadZIP<scalar_type, index_type>::setAbsoluteTolerance(RealT)
    {
      return 0;
    }

    /**
     * @brief Residual contribution of the load is pushed to the bus.
     *
     */
    template <typename scalar_type, typename index_type>
    int LoadZIP<scalar_type, index_type>::evaluateResidual()
    {
      const ScalarT Vr = this->Vr();
      const ScalarT Vi = this->Vi();

      if (alphaI_ == ZERO<RealT> && alphaP_ == ZERO<RealT>)
      {
        Ir() -= G_ * Vr + B_ * Vi;
        Ii() -= G_ * Vi - B_ * Vr;
        return 0;
      }

      const RealT   Vnom2 = Vnom_ * Vnom_;
      const ScalarT V2    = Vr * Vr + Vi * Vi;
      const ScalarT V     = std::sqrt(V2);
      const ScalarT zip   = alphaZ_ + alphaI_ * Vnom_ / V + alphaP_ * Vnom2 / V2;

      Ir() -= (G_ * Vr + B_ * Vi) * zip;
      Ii() -= (G_ * Vi - B_ * Vr) * zip;

      return 0;
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
                      const ScalarT Vr    = this->Vr();
                      const ScalarT Vi    = this->Vi();
                      const RealT   Vnom2 = Vnom_ * Vnom_;
                      const ScalarT V2    = Vr * Vr + Vi * Vi;
                      const ScalarT V     = std::sqrt(V2);
                      const ScalarT zip   = alphaZ_ + alphaI_ * Vnom_ / V + alphaP_ * Vnom2 / V2;
                      return -(G_ * Vr + B_ * Vi) * zip; });
      monitor_->set(Variable::ii, [this]
                    {
                      const ScalarT Vr    = this->Vr();
                      const ScalarT Vi    = this->Vi();
                      const RealT   Vnom2 = Vnom_ * Vnom_;
                      const ScalarT V2    = Vr * Vr + Vi * Vi;
                      const ScalarT V     = std::sqrt(V2);
                      const ScalarT zip   = alphaZ_ + alphaI_ * Vnom_ / V + alphaP_ * Vnom2 / V2;
                      return -(G_ * Vi - B_ * Vr) * zip; });
      monitor_->set(Variable::im, [this]
                    {
                      const ScalarT Vr    = this->Vr();
                      const ScalarT Vi    = this->Vi();
                      const RealT   Vnom2 = Vnom_ * Vnom_;
                      const ScalarT V2    = Vr * Vr + Vi * Vi;
                      const ScalarT V     = std::sqrt(V2);
                      const ScalarT zip   = alphaZ_ + alphaI_ * Vnom_ / V + alphaP_ * Vnom2 / V2;
                      const ScalarT Ir    = -(G_ * Vr + B_ * Vi) * zip;
                      const ScalarT Ii    = -(G_ * Vi - B_ * Vr) * zip;
                      return std::sqrt(Ir * Ir + Ii * Ii); });
      monitor_->set(Variable::p, [this]
                    {
                      const ScalarT Vr    = this->Vr();
                      const ScalarT Vi    = this->Vi();
                      const RealT   Vnom2 = Vnom_ * Vnom_;
                      const ScalarT V2    = Vr * Vr + Vi * Vi;
                      const ScalarT V     = std::sqrt(V2);
                      const ScalarT zip   = alphaZ_ + alphaI_ * Vnom_ / V + alphaP_ * Vnom2 / V2;
                      return -G_ * V2 * zip; });
      monitor_->set(Variable::q, [this]
                    {
                      const ScalarT Vr    = this->Vr();
                      const ScalarT Vi    = this->Vi();
                      const RealT   Vnom2 = Vnom_ * Vnom_;
                      const ScalarT V2    = Vr * Vr + Vi * Vi;
                      const ScalarT V     = std::sqrt(V2);
                      const ScalarT zip   = alphaZ_ + alphaI_ * Vnom_ / V + alphaP_ * Vnom2 / V2;
                      return -B_ * V2 * zip; });
    }

  } // namespace PhasorDynamics
} // namespace GridKit
