#pragma once

#include <GridKit/Model/PhasorDynamics/SignalDelay/SignalDelay.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit::PhasorDynamics
{
  template <class ScalarT, class IdxT>
  SignalDelay<ScalarT, IdxT>::~SignalDelay() = default;

  template <class ScalarT, class IdxT>
  const Model::VariableMonitorBase* SignalDelay<ScalarT, IdxT>::getMonitor() const
  {
    return monitor_.get();
  }

  template <class ScalarT, class IdxT>
  SignalDelay<ScalarT, IdxT>::SignalDelay(const ModelDataT& data, SignalHistory<RealT>& history)
    : delay_(std::get<RealT>(data.parameters.at(SignalDelayParameters::tau))), history_(history)
  {
    if (!std::isfinite(delay_) || delay_ <= 0)
    {
      throw std::invalid_argument("SignalDelay requires a finite positive delay");
    }
    this->size_    = 1;
    this->monitor_ = std::make_unique<MonitorT>(data);
    this->monitor_->set(SignalDelayMonitors::z, [this]
                        { return this->y_.getData()[0]; });
  }

  template <class ScalarT, class IdxT>
  int SignalDelay<ScalarT, IdxT>::setGridKitComponentID(IdxT id)
  {
    this->gridkit_component_id_ = id;
    return 0;
  }

  template <class ScalarT, class IdxT>
  void SignalDelay<ScalarT, IdxT>::linkOutput()
  {
    ports_.out.template port<SignalDelayOutputs::z>().link(this->y_.getData(), &this->getVariableIndex(0));
  }

  template <class ScalarT, class IdxT>
  int SignalDelay<ScalarT, IdxT>::allocate()
  {
    if (!this->allocated_)
    {
      this->allocateVectors(1);
    }
    this->variable_indices_.assign(1, 0);
    this->residual_indices_.assign(1, 0);
    this->tag_.assign(1, false);
    linkOutput();
    this->allocated_ = true;
    return 0;
  }

  template <class ScalarT, class IdxT>
  int SignalDelay<ScalarT, IdxT>::verify() const
  {
    if (!ports_.in.template port<SignalDelayInputs::u>().linked() || !ports_.out.template port<SignalDelayOutputs::z>().linked())
    {
      throw std::invalid_argument("SignalDelay requires linked input and output signals");
    }
    return 0;
  }

  template <class ScalarT, class IdxT>
  int SignalDelay<ScalarT, IdxT>::initialize()
  {
    const auto a           = history_.sample(this->time_, delay_);
    this->y_.getData()[0]  = a.constant + a.weight * ports_.in.template port<SignalDelayInputs::u>().readSignal();
    this->yp_.getData()[0] = 0;
    linkOutput();
    if constexpr (std::is_same_v<ScalarT, DependencyTracking::Variable>)
    {
      this->initializeDependencyTrackingVariableNumbers();
    }
    return 0;
  }

  template <class ScalarT, class IdxT>
  int SignalDelay<ScalarT, IdxT>::tagDifferentiable()
  {
    this->tag_[0] = false;
    return 0;
  }

  template <class ScalarT, class IdxT>
  int SignalDelay<ScalarT, IdxT>::setAbsoluteTolerance(RealT tolerance)
  {
    this->abs_tol_.getData()[0] = tolerance;
    return 0;
  }

  template <class ScalarT, class IdxT>
  int SignalDelay<ScalarT, IdxT>::evaluateResidual()
  {
    const auto a          = history_.sample(this->time_, delay_);
    this->f_.getData()[0] = -this->y_.getData()[0] + a.constant
                            + a.weight * ports_.in.template port<SignalDelayInputs::u>().readSignal();
    return 0;
  }

  template <class ScalarT, class IdxT>
  int SignalDelay<ScalarT, IdxT>::evaluateJacobian()
  {
    if (this->J_rows_buffer_ == nullptr)
    {
      this->nnz_           = 2;
      this->J_rows_buffer_ = new IdxT[2];
      this->J_cols_buffer_ = new IdxT[2];
      this->J_vals_buffer_ = new RealT[2];
    }
    const auto a            = history_.sample(this->time_, delay_);
    const auto source       = ports_.in.template port<SignalDelayInputs::u>().signalVariableIndex();
    this->J_rows_buffer_[0] = this->J_rows_buffer_[1] = this->getResidualIndex(0);
    this->J_cols_buffer_[0]                           = this->getVariableIndex(0);
    this->J_vals_buffer_[0]                           = -1;
    // Preserve the structural entry when the query moves into accepted history.
    this->J_cols_buffer_[1]                           = source == INVALID_INDEX<IdxT> ? this->getVariableIndex(0) : source;
    this->J_vals_buffer_[1]                           = source == INVALID_INDEX<IdxT> ? 0 : a.weight;
    return this->constructCoo();
  }
} // namespace GridKit::PhasorDynamics
