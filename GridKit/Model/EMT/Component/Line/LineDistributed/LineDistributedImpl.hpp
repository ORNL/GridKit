#pragma once

#include <GridKit/Model/EMT/Component/Line/LineDistributed/LineDistributed.hpp>
#include <GridKit/Model/EMT/ComponentInitialization.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    LineDistributed<scalar_type, index_type>::LineDistributed()
    {
      equation_size_ = size_ = static_cast<IdxT>(LineDistributedInternalVariables::MAXIMUM);
      initializePorts();
    }

    template <typename scalar_type, typename index_type>
    LineDistributed<scalar_type, index_type>::LineDistributed(const ModelDataT& data)
      : monitor_(std::make_unique<MonitorT>(data))
    {
      initializeParameters(data);
      equation_size_ = size_ = static_cast<IdxT>(LineDistributedInternalVariables::MAXIMUM);
      propagation_[0]        = std::make_unique<PropagationT>(data.H);
      propagation_[1]        = std::make_unique<PropagationT>(data.H);
      setDerivedParams();
      initializePorts();
      initializeMonitor();
    }

    template <typename scalar_type, typename index_type>
    LineDistributed<scalar_type, index_type>::~LineDistributed()
    {
    }

    template <typename scalar_type, typename index_type>
    void LineDistributed<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Parameter = typename ModelDataT::Parameters;
      N_              = parameter<IdxT>(data, Parameter::N);
      K_              = parameter<IdxT>(data, Parameter::K);
      conductors_     = parameter<ABCVector<IdxT>>(data, Parameter::conductors);
      if (N_ != 3 || K_ != 3 || conductors_ != ABCVector<IdxT>{1, 2, 3})
        throw std::invalid_argument("LineDistributed: expected N=K=3 and conductors=[1,2,3]");
      if (data.Yc.rows != K_ || data.Yc.cols != K_ || data.Yc.validate() || data.H.K != K_)
        throw std::invalid_argument("LineDistributed: invalid coefficient dimensions");
      for (auto pole : data.Yc.poles)
        if (pole.real() >= ZERO<RealT>)
          throw std::invalid_argument("LineDistributed: characteristic admittance must be stable");
      for (const auto& row : data.Yc.E)
        for (auto value : row)
          if (value != ZERO<RealT>)
            throw std::invalid_argument("LineDistributed: characteristic admittance must be proper");
    }

    template <typename scalar_type, typename index_type>
    void LineDistributed<scalar_type, index_type>::setPrehistory(
        RealT omega, const std::array<ABCVector<RealT>, 2>& value, const std::array<ABCVector<RealT>, 2>& derivative)
    {
      if (!std::isfinite(omega) || omega < ZERO<RealT>)
        throw std::invalid_argument("LineDistributed: prehistory frequency must be finite and nonnegative");
      for (size_t e = 0; e < 2; ++e)
        for (size_t p = 0; p < 3; ++p)
          if (!std::isfinite(value[e][p]) || !std::isfinite(derivative[e][p])
              || (omega == ZERO<RealT> && derivative[e][p] != ZERO<RealT>) )
            throw std::invalid_argument("LineDistributed: invalid prehistory values");
      omega_              = omega;
      history_            = value;
      history_derivative_ = derivative;
      has_history_        = true;
    }

    template <typename scalar_type, typename index_type>
    void LineDistributed<scalar_type, index_type>::setDerivedParams()
    {
      for (auto& propagation : propagation_)
      {
        this->addOperator(propagation.get());
        size_ += propagation->size();
      }
    }

    template <typename scalar_type, typename index_type>
    void LineDistributed<scalar_type, index_type>::initializePorts()
    {
      for (size_t e = 0; e < 2; ++e)
      {
        std::vector<SignalT*> input;
        for (size_t p = 0; p < 3; ++p)
        {
          reflected_[e][p].claimProducer();
          input.push_back(&reflected_[e][p]);
          if (propagation_[1 - e])
            signals_.attachSignal(static_cast<LineDistributedExternalVariables>(6 + 3 * e + p), &incidentSignal(e, p));
        }
        if (propagation_[e])
          propagation_[e]->attachInput(input);
      }
    }

    template <typename scalar_type, typename index_type>
    void LineDistributed<scalar_type, index_type>::attachTerminal(size_t end, PhaseSignals characteristic)
    {
      if (allocated_)
        throw std::logic_error("LineDistributed terminals cannot change after allocation");
      if (end > 1)
        throw std::out_of_range("Invalid LineDistributed terminal");
      for (size_t p = 0; p < 3; ++p)
        signals_.attachSignal(static_cast<LineDistributedExternalVariables>(3 * end + p), characteristic[p]);
    }

    template <typename scalar_type, typename index_type>
    typename LineDistributed<scalar_type, index_type>::SignalT& LineDistributed<scalar_type, index_type>::incidentSignal(size_t end, size_t phase)
    {
      if (end > 1 || phase > 2)
        throw std::out_of_range("Invalid LineDistributed incident-current port");
      return propagation_[1 - end]->outputSignal(static_cast<IdxT>(phase));
    }

    template <typename scalar_type, typename index_type>
    typename LineDistributed<scalar_type, index_type>::SignalT& LineDistributed<scalar_type, index_type>::outputSignal(Outputs output)
    {
      const size_t n = static_cast<size_t>(output);
      if (n < 6)
        return reflected_[n / 3][n % 3];
      if (n < 12)
        return incidentSignal((n - 6) / 3, n % 3);
      throw std::out_of_range("LineDistributed: invalid output");
    }

    template <typename scalar_type, typename index_type>
    void LineDistributed<scalar_type, index_type>::assignOutput(Outputs output, SignalT* signal)
    {
      if (allocated_)
        throw std::logic_error("LineDistributed outputs cannot change after allocation");
      if (static_cast<size_t>(output) < static_cast<size_t>(LineDistributedInternalVariables::MAXIMUM))
        signals_.assignSignal(static_cast<LineDistributedInternalVariables>(output), signal);
      else
      {
        auto* value = &outputSignal(output);
        signal->claimProducer();
        signal->setComputed([value]
                            { return value->read(); },
                            [value](typename SignalT::GradientT& gradient, RealT scale)
                            { value->appendGradient(gradient, scale); });
      }
    }

    template <typename scalar_type, typename index_type>
    int LineDistributed<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int LineDistributed<scalar_type, index_type>::allocate()
    {
      if (!allocated_)
        this->allocateVectors(size_);
      const auto size = static_cast<size_t>(size_);
      tag_.resize(size);
      variable_indices_.resize(size);
      residual_indices_.resize(size);
      for (IdxT k = 0; k < static_cast<IdxT>(LineDistributedInternalVariables::MAXIMUM); ++k)
        this->bindSignal(reflected_[static_cast<size_t>(k) / 3][static_cast<size_t>(k) % 3], k);
      const int status = this->allocateOperators();
      if (status != 0)
        return status;
      this->assignGlobalIndices(0);
      this->allocateExternalVectors(static_cast<IdxT>(LineDistributedExternalVariables::MAXIMUM), 0);
      signals_.registerExternalVariableSignals(*this);
      signals_.bindInternalVariableSignals(*this);
      allocated_ = true;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int LineDistributed<scalar_type, index_type>::verify() const
    {
      int error_count = 0;
      if (!signals_.template isAttached<LineDistributedExternalVariables::IC1A>()
          || !signals_.template isAttached<LineDistributedExternalVariables::IC1B>()
          || !signals_.template isAttached<LineDistributedExternalVariables::IC1C>()
          || !signals_.template isAttached<LineDistributedExternalVariables::IC2A>()
          || !signals_.template isAttached<LineDistributedExternalVariables::IC2B>()
          || !signals_.template isAttached<LineDistributedExternalVariables::IC2C>())
      {
        Log::error() << "LineDistributed: a characteristic-admittance input is not attached\n";
        ++error_count;
      }
      for (const auto& propagation : propagation_)
        error_count += propagation ? propagation->verify() : 1;
      return error_count;
    }

    template <typename scalar_type, typename index_type>
    int LineDistributed<scalar_type, index_type>::initialize(const std::map<Outputs, RealT>& outputs)
    {
      this->validateOutputValues(outputs);
      if (!has_history_)
        throw std::invalid_argument("LineDistributed: reflected-current prehistory is required");
      for (size_t e = 0; e < 2; ++e)
      {
        for (size_t p = 0; p < 3; ++p)
        {
          this->checkOutputValue(outputs, static_cast<Outputs>(3 * e + p), history_[e][p]);
          reflected_[e][p].init(static_cast<ScalarT>(history_[e][p]));
          reflected_[e][p].initDerivative(static_cast<ScalarT>(history_derivative_[e][p]));
        }
        const int status = propagation_[e]->initializeSteadyState(omega_, history_[e], history_derivative_[e]);
        if (status != 0)
          return status;
      }
      this->y_.setDataUpdated();
      this->yp_.setDataUpdated();
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int LineDistributed<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      this->abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return this->setAbsoluteToleranceOperators(rel_tol);
    }

    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int LineDistributed<scalar_type, index_type>::evaluateInternalResidual(
        const ScalarT*                  y,
        [[maybe_unused]] const ScalarT* yp,
        const ScalarT*                  y_ext,
        [[maybe_unused]] const ScalarT* yp_ext,
        ScalarT*                        f)
    {
      /* Read variables */
      const ScalarT i_ref1a = y[0];
      const ScalarT i_ref1b = y[1];
      const ScalarT i_ref1c = y[2];
      const ScalarT i_ref2a = y[3];
      const ScalarT i_ref2b = y[4];
      const ScalarT i_ref2c = y[5];

      /* Read external variables */
      const ScalarT i_c1a   = y_ext[0];
      const ScalarT i_c1b   = y_ext[1];
      const ScalarT i_c1c   = y_ext[2];
      const ScalarT i_c2a   = y_ext[3];
      const ScalarT i_c2b   = y_ext[4];
      const ScalarT i_c2c   = y_ext[5];
      const ScalarT i_inc1a = y_ext[6];
      const ScalarT i_inc1b = y_ext[7];
      const ScalarT i_inc1c = y_ext[8];
      const ScalarT i_inc2a = y_ext[9];
      const ScalarT i_inc2b = y_ext[10];
      const ScalarT i_inc2c = y_ext[11];

      /* Reflected-current equations */
      f[0] = -i_ref1a + TWO<RealT> * i_c1a - i_inc1a;
      f[1] = -i_ref1b + TWO<RealT> * i_c1b - i_inc1b;
      f[2] = -i_ref1c + TWO<RealT> * i_c1c - i_inc1c;
      f[3] = -i_ref2a + TWO<RealT> * i_c2a - i_inc2a;
      f[4] = -i_ref2b + TWO<RealT> * i_c2b - i_inc2b;
      f[5] = -i_ref2c + TWO<RealT> * i_c2c - i_inc2c;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int LineDistributed<scalar_type, index_type>::evaluateInternalResidual()
    {
      this->gatherExternalVariables();
      evaluateInternalResidual(y_.getData(), yp_.getData(), y_ext_.data(), yp_ext_.data(), f_.getData());
      const int status = this->evaluateOperatorInternalResiduals();
      f_.setDataUpdated();
      return status;
    }

    template <typename scalar_type, typename index_type>
    int LineDistributed<scalar_type, index_type>::evaluateResidual()
    {
      const int status = evaluateInternalResidual();
      return status == 0 ? this->evaluateExternalResidual() : status;
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* LineDistributed<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    template <typename scalar_type, typename index_type>
    void LineDistributed<scalar_type, index_type>::initializeMonitor()
    {
      using Variable = typename ModelDataT::MonitorableVariables;
      monitor_->set(Variable::i_c1a, [this]
                    { return signals_.template readExternalVariable<LineDistributedExternalVariables::IC1A>(); });
      monitor_->set(Variable::i_c1b, [this]
                    { return signals_.template readExternalVariable<LineDistributedExternalVariables::IC1B>(); });
      monitor_->set(Variable::i_c1c, [this]
                    { return signals_.template readExternalVariable<LineDistributedExternalVariables::IC1C>(); });
      monitor_->set(Variable::i_c2a, [this]
                    { return signals_.template readExternalVariable<LineDistributedExternalVariables::IC2A>(); });
      monitor_->set(Variable::i_c2b, [this]
                    { return signals_.template readExternalVariable<LineDistributedExternalVariables::IC2B>(); });
      monitor_->set(Variable::i_c2c, [this]
                    { return signals_.template readExternalVariable<LineDistributedExternalVariables::IC2C>(); });
      monitor_->set(Variable::i_inc1a, [this]
                    { return incidentSignal(0, 0).read(); });
      monitor_->set(Variable::i_inc1b, [this]
                    { return incidentSignal(0, 1).read(); });
      monitor_->set(Variable::i_inc1c, [this]
                    { return incidentSignal(0, 2).read(); });
      monitor_->set(Variable::i_inc2a, [this]
                    { return incidentSignal(1, 0).read(); });
      monitor_->set(Variable::i_inc2b, [this]
                    { return incidentSignal(1, 1).read(); });
      monitor_->set(Variable::i_inc2c, [this]
                    { return incidentSignal(1, 2).read(); });
      monitor_->set(Variable::i_ref1a, [this]
                    { return y_.getData()[0]; });
      monitor_->set(Variable::i_ref1b, [this]
                    { return y_.getData()[1]; });
      monitor_->set(Variable::i_ref1c, [this]
                    { return y_.getData()[2]; });
      monitor_->set(Variable::i_ref2a, [this]
                    { return y_.getData()[3]; });
      monitor_->set(Variable::i_ref2b, [this]
                    { return y_.getData()[4]; });
      monitor_->set(Variable::i_ref2c, [this]
                    { return y_.getData()[5]; });
    }

  } // namespace EMT
} // namespace GridKit
