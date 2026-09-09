#pragma once

#include <algorithm>
#include <iostream>

#include <GridKit/Model/EMT/Component/Switch/Switch.hpp>
#include <GridKit/Model/EMT/Component/Switch/SwitchData.hpp>
#include <GridKit/Model/EMT/ComponentInitialization.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Constructor for a three-phase switch
     *
     * System sizes:
     * - Number of equations = 3
     * - Number of independent variables = 3
     */
    template <typename scalar_type, typename index_type>
    Switch<scalar_type, index_type>::Switch()
    {
      size_ = 3;
    }

    template <typename scalar_type, typename index_type>
    Switch<scalar_type, index_type>::Switch(const ModelDataT& data)
      : monitor_(std::make_unique<MonitorT>(data))
    {
      initializeParameters(data);
      size_ = 3;
      initializeMonitor();
    }

    template <typename scalar_type, typename index_type>
    Switch<scalar_type, index_type>::~Switch()
    {
    }

    /**
     * @brief Read model parameters from the data object
     */
    template <typename scalar_type, typename index_type>
    void Switch<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Parameter = typename ModelDataT::Parameters;
      i_scale_        = nominalScale<RealT>(data, Parameter::I, std::sqrt(TWO<RealT>));
      if (data.parameters.contains(Parameter::open))
      {
        setOpen(parameter<bool>(data, Parameter::open));
      }
    }

    /**
     * @brief Set the component ID
     */
    template <typename scalar_type, typename index_type>
    int Switch<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /*!
     * @brief allocate method resizes local storage and registers coupling signals.
     */
    template <typename scalar_type, typename index_type>
    int Switch<scalar_type, index_type>::allocate()
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
      this->allocateExternalVectors(static_cast<IdxT>(SwitchExternalVariables::MAXIMUM), 0);
      signals_.registerExternalVariableSignals(*this);

      for (IdxT p = 0; p < 3; ++p)
        this->bindSignal(current_[static_cast<size_t>(p)], p);
      signals_.bindInternalVariableSignals(*this);
      allocated_ = true;
      return 0;
    }

    /**
     * @brief Check model correctness
     *
     * @return Number of model configuration errors found
     */
    template <typename scalar_type, typename index_type>
    int Switch<scalar_type, index_type>::verify() const
    {
      int error_count = 0;

      if (!signals_.template isAttached<SwitchExternalVariables::V1A>()
          || !signals_.template isAttached<SwitchExternalVariables::V1B>()
          || !signals_.template isAttached<SwitchExternalVariables::V1C>()
          || !signals_.template isAttached<SwitchExternalVariables::V2A>()
          || !signals_.template isAttached<SwitchExternalVariables::V2B>()
          || !signals_.template isAttached<SwitchExternalVariables::V2C>())
      {
        Log::error() << "Switch: a terminal voltage port is not attached\n";
        ++error_count;
      }

      if (open_ != ZERO<RealT> && open_ != ONE<RealT>)
      {
        Log::error() << "Switch: the open command must be exactly zero or one\n";
        ++error_count;
      }

      return error_count;
    }

    /**
     * Initialization of the switch model
     *
     * The open command is applied before enforcing the algebraic equations.
     */
    template <typename scalar_type, typename index_type>
    int Switch<scalar_type, index_type>::initialize(const std::map<Outputs, RealT>& outputs)
    {
      this->validateOutputValues(outputs);
      auto* y  = y_.getData();
      auto* yp = yp_.getData();

      y[0]  = this->outputValue(outputs, Outputs::i12a, ZERO<RealT>);
      y[1]  = this->outputValue(outputs, Outputs::i12b, ZERO<RealT>);
      y[2]  = this->outputValue(outputs, Outputs::i12c, ZERO<RealT>);
      yp[0] = 0.0;
      yp[1] = 0.0;
      yp[2] = 0.0;

      y_.setDataUpdated();
      yp_.setDataUpdated();

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
    int Switch<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol * i_scale_));
      return 0;
    }

    /**
     * @brief Internal residual
     *
     * The open command is a constant mask during differentiation, so the
     * open and closed configurations share one residual row structure.
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int Switch<scalar_type, index_type>::evaluateInternalResidual(
        const ScalarT*                  y,
        [[maybe_unused]] const ScalarT* yp,
        const ScalarT*                  y_ext,
        [[maybe_unused]] const ScalarT* yp_ext,
        ScalarT*                        f)
    {
      /* Read variables */
      const ScalarT i12a = y[0];
      const ScalarT i12b = y[1];
      const ScalarT i12c = y[2];

      // Set coupling variable aliases
      const ScalarT v1a = y_ext[0];
      const ScalarT v1b = y_ext[1];
      const ScalarT v1c = y_ext[2];
      const ScalarT v2a = y_ext[3];
      const ScalarT v2b = y_ext[4];
      const ScalarT v2c = y_ext[5];

      const RealT closed = ONE<RealT> - open_;

      /* 3 switch algebraic equations */
      f[0] = open_ * i12a + closed * (v2a - v1a);
      f[1] = open_ * i12b + closed * (v2b - v1b);
      f[2] = open_ * i12c + closed * (v2c - v1c);

      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Switch<scalar_type, index_type>::evaluateInternalResidual()
    {
      this->gatherExternalVariables();

      const auto* y  = y_.getData();
      const auto* yp = yp_.getData();
      auto*       f  = f_.getData();
      evaluateInternalResidual(y, yp, y_ext_.data(), yp_ext_.data(), f);
      f_.setDataUpdated();

      return 0;
    }

    /**
     * @brief Assemble the switch equations.
     *
     */
    template <typename scalar_type, typename index_type>
    int Switch<scalar_type, index_type>::evaluateResidual()
    {
      return evaluateInternalResidual();
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* Switch<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    template <typename scalar_type, typename index_type>
    void Switch<scalar_type, index_type>::initializeMonitor()
    {
      using Variable = typename ModelDataT::MonitorableVariables;

      monitor_->set(Variable::open, [this]
                    { return static_cast<ScalarT>(open_); });
      monitor_->set(Variable::i12a, [this]
                    { return y_.getData()[0]; });
      monitor_->set(Variable::i12b, [this]
                    { return y_.getData()[1]; });
      monitor_->set(Variable::i12c, [this]
                    { return y_.getData()[2]; });
    }

  } // namespace EMT
} // namespace GridKit
