#pragma once

#include <cmath>

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
     * System sizes:
     * - Number of equations = 2
     * - Number of independent variables = 2
     */
    template <typename scalar_type, typename index_type>
    LoadZIP<scalar_type, index_type>::LoadZIP(BusT* bus)
      : bus_(bus)
    {
      size_ = 2;
      setDerivedParams();
    }

    template <typename scalar_type, typename index_type>
    LoadZIP<scalar_type, index_type>::LoadZIP(BusT* bus, RealT Pnom, RealT Qnom, RealT Vnom, RealT alphaI, RealT alphaP)
      : bus_(bus),
        Pnom_(Pnom),
        Qnom_(Qnom),
        Vnom_(Vnom),
        alphaI_(alphaI),
        alphaP_(alphaP)
    {
      size_ = 2;
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
      size_ = 2;
    }

    template <typename scalar_type, typename index_type>
    LoadZIP<scalar_type, index_type>::~LoadZIP()
    {
    }

    template <typename scalar_type, typename index_type>
    void LoadZIP<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      if (data.parameters.contains(ModelDataT::Parameters::Pnom))
      {
        Pnom_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::Pnom));
      }

      if (data.parameters.contains(ModelDataT::Parameters::Qnom))
      {
        Qnom_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::Qnom));
      }

      if (data.parameters.contains(ModelDataT::Parameters::Vnom))
      {
        Vnom_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::Vnom));
      }

      if (data.parameters.contains(ModelDataT::Parameters::alphaI))
      {
        alphaI_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::alphaI));
      }

      if (data.parameters.contains(ModelDataT::Parameters::alphaP))
      {
        alphaP_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::alphaP));
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
      auto size = static_cast<size_t>(size_); // avoid compiler warnings
      f_.resize(size);
      y_.resize(size);
      yp_.resize(size);
      tag_.resize(size);
      abs_tol_.resize(size);
      variable_indices_.resize(size);
      residual_indices_.resize(size);

      // Resize coupling data
      wb_.resize(2);
      h_.resize(2);

      // Default variable and residual index mapping to local index
      for (IdxT j = 0; j < size_; ++j)
      {
        this->setVariableIndex(j, j);
        this->setResidualIndex(j, j);
      }

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

      const RealT   Vnom2 = Vnom_ * Vnom_;
      const ScalarT V2    = vr * vr + vi * vi;
      const ScalarT V     = std::sqrt(V2);
      const ScalarT zip   = alphaZ_ + alphaI_ * Vnom_ / V + alphaP_ * Vnom2 / V2;

      y_[0] = -(G_ * vr + B_ * vi) * zip;
      y_[1] = -(G_ * vi - B_ * vr) * zip;

      yp_[0] = 0.0;
      yp_[1] = 0.0;
      return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <typename scalar_type, typename index_type>
    int LoadZIP<scalar_type, index_type>::tagDifferentiable()
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
    int LoadZIP<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      std::fill(abs_tol_.begin(), abs_tol_.end(), rel_tol);
      return 0;
    }

    /**
     * @brief Bus residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int LoadZIP<scalar_type, index_type>::evaluateBusResidual(
        const ScalarT*                  y,
        [[maybe_unused]] const ScalarT* yp,
        [[maybe_unused]] const ScalarT* wb,
        ScalarT*                        h)
    {
      const ScalarT Ir = y[0];
      const ScalarT Ii = y[1];
      h[0]             = Ir;
      h[1]             = Ii;

      return 0;
    }

    /**
     * @brief Residual contribution of the load is pushed to the bus.
     *
     */
    template <typename scalar_type, typename index_type>
    int LoadZIP<scalar_type, index_type>::evaluateResidual()
    {
      wb_[0] = Vr();
      wb_[1] = Vi();
      evaluateInternalResidual(y_.data(), yp_.data(), wb_.data(), f_.data());
      evaluateBusResidual(y_.data(), yp_.data(), wb_.data(), h_.data());
      Ir() += h_[0];
      Ii() += h_[1];

      return 0;
    }

    /**
     * @brief Internal residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int LoadZIP<scalar_type, index_type>::evaluateInternalResidual(
        const ScalarT*                  y,
        [[maybe_unused]] const ScalarT* yp,
        const ScalarT*                  wb,
        ScalarT*                        f)
    {
      const ScalarT Vr    = wb[0];
      const ScalarT Vi    = wb[1];
      const ScalarT Ir    = y[0];
      const ScalarT Ii    = y[1];
      const RealT   Vnom2 = Vnom_ * Vnom_;
      const ScalarT V2    = Vr * Vr + Vi * Vi;
      const ScalarT V     = std::sqrt(V2);
      const ScalarT zip   = alphaZ_ + alphaI_ * Vnom_ / V + alphaP_ * Vnom2 / V2;

      f[0] = Ir + (G_ * Vr + B_ * Vi) * zip;
      f[1] = Ii + (G_ * Vi - B_ * Vr) * zip;

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
                    { return y_[0]; });
      monitor_->set(Variable::ii, [this]
                    { return y_[1]; });
      monitor_->set(Variable::im, [this]
                    { return std::sqrt(y_[0] * y_[0] + y_[1] * y_[1]); });
      monitor_->set(Variable::p, [this]
                    { return Vr() * y_[0] + Vi() * y_[1]; });
      monitor_->set(Variable::q, [this]
                    { return Vi() * y_[0] - Vr() * y_[1]; });
    }

  } // namespace PhasorDynamics
} // namespace GridKit
