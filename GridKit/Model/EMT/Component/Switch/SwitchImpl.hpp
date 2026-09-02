#pragma once

#include <algorithm>
#include <iostream>

#include <GridKit/Model/EMT/Component/Switch/Switch.hpp>
#include <GridKit/Model/EMT/Component/Switch/SwitchData.hpp>
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
      if (data.parameters.contains(Parameter::open))
      {
        setOpen(std::get<bool>(data.parameters.at(Parameter::open)));
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
     * @brief allocate method resizes local storage and registers coupling nodes.
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
      this->allocateExternalVectors(static_cast<IdxT>(SwitchExternalVariables::MAXIMUM), 6);
      signals_.registerExternalVariableNodes(*this);
      this->setExternalResidualNode(0, signals_.template getAttachedSignalNode<SwitchExternalVariables::V1A>());
      this->setExternalResidualNode(1, signals_.template getAttachedSignalNode<SwitchExternalVariables::V1B>());
      this->setExternalResidualNode(2, signals_.template getAttachedSignalNode<SwitchExternalVariables::V1C>());
      this->setExternalResidualNode(3, signals_.template getAttachedSignalNode<SwitchExternalVariables::V2A>());
      this->setExternalResidualNode(4, signals_.template getAttachedSignalNode<SwitchExternalVariables::V2B>());
      this->setExternalResidualNode(5, signals_.template getAttachedSignalNode<SwitchExternalVariables::V2C>());

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
    int Switch<scalar_type, index_type>::initialize()
    {
      auto* y  = y_.getData();
      auto* yp = yp_.getData();

      y[0]  = 0.0;
      y[1]  = 0.0;
      y[2]  = 0.0;
      yp[0] = 0.0;
      yp[1] = 0.0;
      yp[2] = 0.0;

      y_.setDataUpdated();
      yp_.setDataUpdated();

      return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <typename scalar_type, typename index_type>
    int Switch<scalar_type, index_type>::tagDifferentiable()
    {
      tag_[0] = false;
      tag_[1] = false;
      tag_[2] = false;

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
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
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

    /**
     * @brief External residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int Switch<scalar_type, index_type>::evaluateExternalResidual(
        const ScalarT*                  y,
        [[maybe_unused]] const ScalarT* yp,
        [[maybe_unused]] const ScalarT* y_ext,
        [[maybe_unused]] const ScalarT* yp_ext,
        ScalarT*                        f_ext)
    {
      const ScalarT i12a = y[0];
      const ScalarT i12b = y[1];
      const ScalarT i12c = y[2];

      f_ext[0] = -i12a;
      f_ext[1] = -i12b;
      f_ext[2] = -i12c;
      f_ext[3] = i12a;
      f_ext[4] = i12b;
      f_ext[5] = i12c;

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
     * @brief External residual contributions to the terminal buses.
     *
     */
    template <typename scalar_type, typename index_type>
    int Switch<scalar_type, index_type>::evaluateExternalResidual()
    {
      const auto* y  = y_.getData();
      const auto* yp = yp_.getData();
      evaluateExternalResidual(y, yp, y_ext_.data(), yp_ext_.data(), f_ext_.data());
      this->scatterExternalResidual();

      return 0;
    }

    /**
     * @brief Residual contribution of the switch is pushed to the buses.
     *
     */
    template <typename scalar_type, typename index_type>
    int Switch<scalar_type, index_type>::evaluateResidual()
    {
      evaluateInternalResidual();
      return evaluateExternalResidual();
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
