#pragma once

#include <algorithm>
#include <iostream>

#include <GridKit/Model/PhasorDynamics/Load/LoadZ/LoadZ.hpp>
#include <GridKit/Model/PhasorDynamics/Load/LoadZ/LoadZData.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Constructor for a constant-impedance load
     *
     * System sizes:
     * - Number of equations = 2
     * - Number of independent variables = 2
     */
    template <typename scalar_type, typename index_type>
    LoadZ<scalar_type, index_type>::LoadZ()
    {
      size_ = 2;
      setDerivedParams();
    }

    template <typename scalar_type, typename index_type>
    LoadZ<scalar_type, index_type>::LoadZ(RealT R, RealT X)
      : R_(R),
        X_(X)
    {
      size_ = 2;
      setDerivedParams();
    }

    template <typename scalar_type, typename index_type>
    LoadZ<scalar_type, index_type>::LoadZ(const ModelDataT& data)
      : monitor_(std::make_unique<MonitorT>(data))
    {
      using Parameter = typename ModelDataT::Parameters;
      if (data.parameters.contains(Parameter::R))
      {
        R_ = std::get<RealT>(data.parameters.at(Parameter::R));
      }

      if (data.parameters.contains(Parameter::X))
      {
        X_ = std::get<RealT>(data.parameters.at(Parameter::X));
      }

      size_ = 2;
      setDerivedParams();
      initializeMonitor();
    }

    template <typename scalar_type, typename index_type>
    LoadZ<scalar_type, index_type>::~LoadZ()
    {
    }

    /**
     * @brief Set the component ID
     */
    template <typename scalar_type, typename index_type>
    int LoadZ<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /*!
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <typename scalar_type, typename index_type>
    int LoadZ<scalar_type, index_type>::allocate()
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
      this->allocateExternalVectors(static_cast<IdxT>(Utilities::enum_size<LoadZExternalVariables>()));

      allocated_ = true;
      return 0;
    }

    /**
     * Initialization of the load model
     *
     */
    template <typename scalar_type, typename index_type>
    int LoadZ<scalar_type, index_type>::initialize()
    {
      ScalarT vr = Vr();
      ScalarT vi = Vi();
      ScalarT ir = -(g_ * vr - b_ * vi);
      ScalarT ii = -(b_ * vr + g_ * vi);

      auto* y  = y_.getData();
      auto* yp = yp_.getData();

      y[0]  = ir;
      y[1]  = ii;
      yp[0] = 0.0;
      yp[1] = 0.0;

      y_.setDataUpdated();
      yp_.setDataUpdated();

      // For DependencyTracking::Variable, set variable numbers
      if constexpr (std::is_same_v<scalar_type, DependencyTracking::Variable>)
      {
        this->initializeDependencyTrackingVariableNumbers();
      }

      return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <typename scalar_type, typename index_type>
    int LoadZ<scalar_type, index_type>::tagDifferentiable()
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
    int LoadZ<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return 0;
    }

    /**
     * @brief External residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int LoadZ<scalar_type, index_type>::evaluateExternalResidual(
        const ScalarT*                  y,
        [[maybe_unused]] const ScalarT* yp,
        [[maybe_unused]] const ScalarT* y_ext,
        ScalarT*                        f_ext)
    {
      const auto VR = static_cast<size_t>(LoadZExternalVariables::VR);
      const auto VI = static_cast<size_t>(LoadZExternalVariables::VI);

      const ScalarT Ir = y[0];
      const ScalarT Ii = y[1];
      f_ext[VR]        = Ir;
      f_ext[VI]        = Ii;

      return 0;
    }

    /**
     * @brief Internal residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int LoadZ<scalar_type, index_type>::evaluateInternalResidual(
        const ScalarT*                  y,
        [[maybe_unused]] const ScalarT* yp,
        const ScalarT*                  y_ext,
        ScalarT*                        f)
    {
      const auto VR = static_cast<size_t>(LoadZExternalVariables::VR);
      const auto VI = static_cast<size_t>(LoadZExternalVariables::VI);

      const ScalarT Vr = y_ext[VR];
      const ScalarT Vi = y_ext[VI];
      const ScalarT Ir = y[0];
      const ScalarT Ii = y[1];
      f[0]             = Ir + g_ * Vr - b_ * Vi;
      f[1]             = Ii + b_ * Vr + g_ * Vi;

      return 0;
    }

    /**
     * @brief Gather external variables and index maps.
     *
     */
    template <typename scalar_type, typename index_type>
    void LoadZ<scalar_type, index_type>::gatherExternalVariables()
    {
      const auto VR = static_cast<size_t>(LoadZExternalVariables::VR);
      const auto VI = static_cast<size_t>(LoadZExternalVariables::VI);

      auto* y_ext   = y_ext_.getData();
      auto* columns = variable_indices_ext_.data();
      auto* rows    = residual_indices_ext_.data();

      ports_.terminal_vr[LoadZBuses::bus].gather(y_ext[VR], columns[VR], rows[VR]);
      ports_.terminal_vi[LoadZBuses::bus].gather(y_ext[VI], columns[VI], rows[VI]);
    }

    /**
     * @brief Evaluate the residual rows the load owns.
     */
    template <typename scalar_type, typename index_type>
    int LoadZ<scalar_type, index_type>::evaluateInternalResidual()
    {
      return evaluateInternalResidual(y_.getData(), yp_.getData(), y_ext_.getData(), f_.getData());
    }

    /**
     * @brief Evaluate the load current injection into the bus.
     */
    template <typename scalar_type, typename index_type>
    int LoadZ<scalar_type, index_type>::evaluateExternalResidual()
    {
      return evaluateExternalResidual(y_.getData(), yp_.getData(), y_ext_.getData(), f_ext_.getData());
    }

    /**
     * @brief Derived parameters
     *
     */
    template <typename scalar_type, typename index_type>
    void LoadZ<scalar_type, index_type>::setDerivedParams()
    {
      b_ = -X_ / (R_ * R_ + X_ * X_);
      g_ = R_ / (R_ * R_ + X_ * X_);
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* LoadZ<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    template <typename scalar_type, typename index_type>
    void LoadZ<scalar_type, index_type>::initializeMonitor()
    {
      using Variable = typename ModelDataT::MonitorableVariables;

      monitor_->set(Variable::p, [this]
                    { return Vr() * y_.getData()[0] + Vi() * y_.getData()[1]; });
      monitor_->set(Variable::q, [this]
                    { return Vi() * y_.getData()[0] - Vr() * y_.getData()[1]; });
    }

  } // namespace PhasorDynamics
} // namespace GridKit
