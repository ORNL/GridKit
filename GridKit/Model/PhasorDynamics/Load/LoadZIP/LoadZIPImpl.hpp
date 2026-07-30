#pragma once

#include <cmath>

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
    LoadZIP<scalar_type, index_type>::LoadZIP()
    {
      size_ = 2;
      setDerivedParams();
    }

    template <typename scalar_type, typename index_type>
    LoadZIP<scalar_type, index_type>::LoadZIP(RealT Pnom, RealT Qnom, RealT Vnom, RealT alphaI, RealT alphaP)
      : Pnom_(Pnom),
        Qnom_(Qnom),
        Vnom_(Vnom),
        alphaI_(alphaI),
        alphaP_(alphaP)
    {
      size_ = 2;
      setDerivedParams();
    }

    template <typename scalar_type, typename index_type>
    LoadZIP<scalar_type, index_type>::LoadZIP(const ModelDataT& data)
      : monitor_(std::make_unique<MonitorT>(data))
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
      using Parameter = typename ModelDataT::Parameters;
      if (data.parameters.contains(Parameter::Pnom))
      {
        Pnom_ = std::get<RealT>(data.parameters.at(Parameter::Pnom));
      }

      if (data.parameters.contains(Parameter::Qnom))
      {
        Qnom_ = std::get<RealT>(data.parameters.at(Parameter::Qnom));
      }

      if (data.parameters.contains(Parameter::Vnom))
      {
        Vnom_ = std::get<RealT>(data.parameters.at(Parameter::Vnom));
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
      this->allocateExternalVectors(static_cast<IdxT>(LoadZIPExternalVariables::MAXIMUM));
      f_ext_.resize(2);
      residual_indices_ext_.assign(2, INVALID_INDEX<IdxT>);

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

      const RealT   Vnom2 = Vnom_ * Vnom_;
      const ScalarT V2    = vr * vr + vi * vi;
      const ScalarT V     = std::sqrt(V2);
      const ScalarT zip   = alphaZ_ + alphaI_ * Vnom_ / V + alphaP_ * Vnom2 / V2;

      auto* y  = y_.getData();
      auto* yp = yp_.getData();

      y[0]  = -(G_ * vr + B_ * vi) * zip;
      y[1]  = -(G_ * vi - B_ * vr) * zip;
      yp[0] = 0.0;
      yp[1] = 0.0;
      y_.setDataUpdated();
      yp_.setDataUpdated();
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
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return 0;
    }

    /**
     * @brief External residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int LoadZIP<scalar_type, index_type>::evaluateExternalResidual(
        const ScalarT*                  y,
        [[maybe_unused]] const ScalarT* yp,
        [[maybe_unused]] const ScalarT* y_ext,
        ScalarT*                        f_ext)
    {
      const ScalarT Ir = y[0];
      const ScalarT Ii = y[1];
      f_ext[0]         = Ir;
      f_ext[1]         = Ii;

      return 0;
    }

    /**
     * @brief Gather external variables and index maps.
     *
     */
    template <typename scalar_type, typename index_type>
    void LoadZIP<scalar_type, index_type>::gatherExternalVariables()
    {
      static constexpr auto VR = LoadZIPExternalVariables::VR;
      static constexpr auto VI = LoadZIPExternalVariables::VI;

      y_ext_[0]                = Vr();
      y_ext_[1]                = Vi();
      variable_indices_ext_[0] = signals_.template readExternalVariableIndex<VR>();
      variable_indices_ext_[1] = signals_.template readExternalVariableIndex<VI>();
      residual_indices_ext_[0] = signals_.template readExternalResidualIndex<VR>();
      residual_indices_ext_[1] = signals_.template readExternalResidualIndex<VI>();
    }

    /**
     * @brief Internal residual for the load model.
     *
     */
    template <typename scalar_type, typename index_type>
    int LoadZIP<scalar_type, index_type>::evaluateInternalResidual()
    {
      gatherExternalVariables();

      const auto* y  = y_.getData();
      const auto* yp = yp_.getData();
      auto*       f  = f_.getData();
      evaluateInternalResidual(y, yp, y_ext_.data(), f);
      f_.setDataUpdated();

      return 0;
    }

    /**
     * @brief External residual contributions to rows owned by the bus.
     *
     */
    template <typename scalar_type, typename index_type>
    int LoadZIP<scalar_type, index_type>::evaluateExternalResidual()
    {
      const auto* y  = y_.getData();
      const auto* yp = yp_.getData();
      evaluateExternalResidual(y, yp, y_ext_.data(), f_ext_.data());

      return 0;
    }

    /**
     * @brief Evaluate the internal residual and external residual
     * contributions.
     *
     */
    template <typename scalar_type, typename index_type>
    int LoadZIP<scalar_type, index_type>::evaluateResidual()
    {
      evaluateInternalResidual();
      evaluateExternalResidual();

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
        const ScalarT*                  y_ext,
        ScalarT*                        f)
    {
      const ScalarT Vr    = y_ext[0];
      const ScalarT Vi    = y_ext[1];
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
                    { return y_.getData()[0]; });
      monitor_->set(Variable::ii, [this]
                    { return y_.getData()[1]; });
      monitor_->set(Variable::im, [this]
                    { return std::sqrt(y_.getData()[0] * y_.getData()[0] + y_.getData()[1] * y_.getData()[1]); });
      monitor_->set(Variable::p, [this]
                    { return Vr() * y_.getData()[0] + Vi() * y_.getData()[1]; });
      monitor_->set(Variable::q, [this]
                    { return Vi() * y_.getData()[0] - Vr() * y_.getData()[1]; });
    }

  } // namespace PhasorDynamics
} // namespace GridKit
