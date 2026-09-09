#pragma once

#include <cmath>
#include <limits>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Component/Controller/InnerCurrentControl/InnerCurrentControl.hpp>
#include <GridKit/Model/EMT/ComponentInitialization.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template <typename scalar_type, typename index_type>
      InnerCurrentControl<scalar_type, index_type>::InnerCurrentControl()
        : InnerCurrentControl(ModelDataT{})
      {
      }

      template <typename scalar_type, typename index_type>
      InnerCurrentControl<scalar_type, index_type>::InnerCurrentControl(const ModelDataT& data)
        : monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        size_ = 8;
        signals_.template assignSignal<InnerCurrentControlInternalVariables::ILIMD>(&output_[0]);
        signals_.template assignSignal<InnerCurrentControlInternalVariables::ILIMQ>(&output_[1]);
        signals_.template assignSignal<InnerCurrentControlInternalVariables::UD>(&output_[2]);
        signals_.template assignSignal<InnerCurrentControlInternalVariables::UQ>(&output_[3]);
        initializeMonitor();
      }

      template <typename scalar_type, typename index_type>
      InnerCurrentControl<scalar_type, index_type>::~InnerCurrentControl() = default;

      template <typename scalar_type, typename index_type>
      void InnerCurrentControl<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Parameter = typename ModelDataT::Parameters;
        i_scale_        = nominalScale<RealT>(data, Parameter::I, std::sqrt(THREE<RealT>));
        v_scale_        = nominalScale<RealT>(data, Parameter::V, ONE<RealT>);
        C_              = parameter<RealT>(data, Parameter::C, C_);
        Tf_             = parameter<RealT>(data, Parameter::Tf, Tf_);
        L_              = parameter<RealT>(data, Parameter::L, L_);
        Kp_             = parameter<RealT>(data, Parameter::Kp, Kp_);
        Ki_             = parameter<RealT>(data, Parameter::Ki, Ki_);
        Kaw_            = parameter<RealT>(data, Parameter::Kaw, Kaw_);
        Imax_           = parameter<RealT>(data, Parameter::Imax, Imax_);
        if (Imax_ > ZERO<RealT>)
        {
          ai_ = ONE<RealT> / (Imax_ * Imax_);
        }
      }

      template <typename scalar_type, typename index_type>
      void InnerCurrentControl<scalar_type, index_type>::attachInput(InputSignals inputs)
      {
        if (allocated_)
          throw std::logic_error("InnerCurrentControl inputs cannot change after allocation");
        for (size_t p = 0; p < inputs.size(); ++p)
          signals_.attachSignal(static_cast<InnerCurrentControlExternalVariables>(p), inputs[p]);
      }

      template <typename scalar_type, typename index_type>
      void InnerCurrentControl<scalar_type, index_type>::assignOutput(Outputs output, SignalT* signal)
      {
        if (allocated_)
          throw std::logic_error("InnerCurrentControl outputs cannot change after allocation");
        const auto n = static_cast<size_t>(output);
        if (n >= output_.size() || !signal || alias_[n])
          throw std::invalid_argument("InnerCurrentControl: invalid or duplicate output assignment");
        signal->claimProducer();
        alias_[n] = signal;
      }

      template <typename scalar_type, typename index_type>
      typename InnerCurrentControl<scalar_type, index_type>::SignalT& InnerCurrentControl<scalar_type, index_type>::outputSignal(Outputs output)
      {
        return output_.at(static_cast<size_t>(output));
      }

      template <typename scalar_type, typename index_type>
      int InnerCurrentControl<scalar_type, index_type>::setGridKitComponentID(IdxT id)
      {
        gridkit_component_id_ = id;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int InnerCurrentControl<scalar_type, index_type>::allocate()
      {
        if (!allocated_)
          this->allocateVectors(size_);
        tag_.resize(static_cast<size_t>(size_));
        variable_indices_.resize(static_cast<size_t>(size_));
        residual_indices_.resize(static_cast<size_t>(size_));
        this->assignGlobalIndices(0);
        this->allocateExternalVectors(static_cast<IdxT>(InnerCurrentControlExternalVariables::MAXIMUM), 0);
        signals_.registerExternalVariableSignals(*this);
        signals_.bindInternalVariableSignals(*this);
        for (size_t n = 0; n < alias_.size(); ++n)
          if (alias_[n])
            this->bindSignal(*alias_[n], static_cast<IdxT>(2 + n));
        allocated_ = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int InnerCurrentControl<scalar_type, index_type>::verify() const
      {
        int error_count   = 0;
        using V           = InnerCurrentControlExternalVariables;
        const auto inputs = signals_.attachedSignals({V::VD, V::VQ, V::ID, V::IQ, V::ICMDD, V::ICMDQ, V::OMEGA, V::ULIMD, V::ULIMQ});
        if (inputs.size() != static_cast<size_t>(V::MAXIMUM))
        {
          Log::error() << "InnerCurrentControl: all inputs are required\n";
          ++error_count;
        }
        for (const auto* input : inputs)
          if (!input->linked())
          {
            Log::error() << "InnerCurrentControl: inputs must have linked sources\n";
            ++error_count;
          }
        if (!std::isfinite(Tf_) || Tf_ <= ZERO<RealT>)
        {
          Log::error() << "InnerCurrentControl: Tf must be finite and positive\n";
          ++error_count;
        }
        if (!std::isfinite(C_) || C_ < ZERO<RealT>)
        {
          Log::error() << "InnerCurrentControl: C must be finite and nonnegative\n";
          ++error_count;
        }
        for (const auto value : {L_, Kp_, Ki_, Kaw_, Imax_, ai_})
          if (!std::isfinite(value) || value <= ZERO<RealT>)
          {
            Log::error() << "InnerCurrentControl: parameters must be finite and positive\n";
            ++error_count;
          }
        return error_count;
      }

      template <typename scalar_type, typename index_type>
      void InnerCurrentControl<scalar_type, index_type>::validateInitialState(const std::map<std::string, RealT>& values) const
      {
        this->template parseInitialOutputs<InnerCurrentControl>(values);
      }

      template <typename scalar_type, typename index_type>
      int InnerCurrentControl<scalar_type, index_type>::initializeState(const std::map<std::string, RealT>& values)
      {
        return this->initializeOutputs(*this, values);
      }

      template <typename scalar_type, typename index_type>
      int InnerCurrentControl<scalar_type, index_type>::initialize(const std::map<Outputs, RealT>& outputs)
      {
        this->validateOutputValues(outputs);
        if (const int errors = verify(); errors != 0)
          return errors;
        using V           = InnerCurrentControlExternalVariables;
        const RealT vd    = static_cast<RealT>(signals_.template readExternalVariable<V::VD>());
        const RealT vq    = static_cast<RealT>(signals_.template readExternalVariable<V::VQ>());
        const RealT id    = static_cast<RealT>(signals_.template readExternalVariable<V::ID>());
        const RealT iq    = static_cast<RealT>(signals_.template readExternalVariable<V::IQ>());
        const RealT omega = static_cast<RealT>(signals_.template readExternalVariable<V::OMEGA>());
        const RealT icmdd = static_cast<RealT>(signals_.template readExternalVariable<V::ICMDD>()) - omega * C_ * vq;
        const RealT icmdq = static_cast<RealT>(signals_.template readExternalVariable<V::ICMDQ>()) + omega * C_ * vd;
        for (const auto value : {vd, vq, id, iq, icmdd, icmdq, omega})
          if (!std::isfinite(value))
            throw std::invalid_argument("InnerCurrentControl: nonfinite initial input");
        const RealT                li = std::sqrt(Math::max(ONE<RealT>, ai_ * (icmdd * icmdd + icmdq * icmdq)));
        const std::array<RealT, 2> limited{icmdd / li, icmdq / li};
        const RealT                tolerance = RealT{64} * std::numeric_limits<RealT>::epsilon();
        for (size_t n = 0; n < 2; ++n)
        {
          const auto  key   = static_cast<Outputs>(n);
          const RealT value = this->outputValue(outputs, key, limited[n]);
          if (std::abs(value - limited[n]) > tolerance * (ONE<RealT> + std::abs(limited[n])))
            throw std::invalid_argument("InnerCurrentControl: initial limited current must match the current command");
        }
        const std::array<RealT, 2> base{vd - omega * L_ * iq + Kp_ * (limited[0] - id),
                                        vq + omega * L_ * id + Kp_ * (limited[1] - iq)};
        const std::array<RealT, 2> u{this->outputValue(outputs, Outputs::ud, base[0]),
                                     this->outputValue(outputs, Outputs::uq, base[1])};
        const std::array<RealT, 6> values{u[0] - base[0], u[1] - base[1], limited[0], limited[1], u[0], u[1]};
        for (const RealT value : values)
          if (!std::isfinite(value))
            throw std::invalid_argument("InnerCurrentControl: nonfinite derived initial value");
        auto* y = y_.getData();
        for (size_t n = 0; n < values.size(); ++n)
          y[n] = values[n];
        y[6]     = vd;
        y[7]     = vq;
        auto* yp = yp_.getData();
        for (IdxT n = 0; n < size_; ++n)
          yp[n] = ZERO<RealT>;
        y_.setDataUpdated();
        yp_.setDataUpdated();
        return 0;
      }

      template <typename scalar_type, typename index_type>
      typename Component<scalar_type, index_type>::InitializationPortsT InnerCurrentControl<scalar_type, index_type>::initializationPorts()
      {
        using V = InnerCurrentControlExternalVariables;
        typename Component<ScalarT, IdxT>::InitializationPortsT ports;
        ports.inputs  = signals_.attachedSignals({V::VD, V::VQ, V::ID, V::IQ, V::ICMDD, V::ICMDQ, V::OMEGA});
        ports.targets = signals_.attachedSignals({V::ICMDD, V::ICMDQ});
        for (size_t n = 0; n < output_.size(); ++n)
        {
          const auto name = std::string(magic_enum::enum_name(static_cast<Outputs>(n)));
          ports.outputs.emplace(name, &output_[n]);
          if (alias_[n])
            ports.outputs.emplace(name, alias_[n]);
        }
        return ports;
      }

      template <typename scalar_type, typename index_type>
      void InnerCurrentControl<scalar_type, index_type>::prepareInitialization(typename Component<ScalarT, IdxT>::InitialStateT& initial)
      {
        if (initial.omega() == ZERO<RealT>)
          return;
        using V        = InnerCurrentControlExternalVariables;
        const RealT id = initial.value(*signals_.template getAttachedSignal<V::ID>());
        const RealT iq = initial.value(*signals_.template getAttachedSignal<V::IQ>());
        const RealT li = std::sqrt(Math::max(ONE<RealT>, ai_ * (id * id + iq * iq)));
        if (std::abs(li - ONE<RealT>) > RealT{1e-10})
          throw std::invalid_argument("InnerCurrentControl: initial current must lie inside the current limit");
        const RealT vd = initial.value(*signals_.template getAttachedSignal<V::VD>());
        const RealT vq = initial.value(*signals_.template getAttachedSignal<V::VQ>());
        initial.require(*signals_.template getAttachedSignal<V::ICMDD>(), id + initial.omega() * C_ * vq, *this);
        initial.require(*signals_.template getAttachedSignal<V::ICMDQ>(), iq - initial.omega() * C_ * vd, *this);
        initial.provide(output_[0], id / li);
        initial.provide(output_[1], iq / li);
        const auto outputs = this->template parseInitialOutputs<InnerCurrentControl>(initial.outputs(*this));
        for (size_t n = 2; n < output_.size(); ++n)
        {
          const auto key = static_cast<Outputs>(n);
          if (!outputs.contains(key))
            throw std::invalid_argument("InnerCurrentControl: balanced initialization requires ud and uq");
          initial.provide(output_[n], outputs.at(key));
        }
      }

      template <typename scalar_type, typename index_type>
      int InnerCurrentControl<scalar_type, index_type>::setAbsoluteTolerance(RealT tolerance)
      {
        auto* absolute = abs_tol_.getData();
        for (size_t p = 0; p < 2; ++p)
        {
          absolute[p]     = tolerance * v_scale_;
          absolute[2 + p] = tolerance * i_scale_;
          absolute[4 + p] = tolerance * v_scale_;
          absolute[6 + p] = tolerance * v_scale_;
        }
        abs_tol_.setDataUpdated();
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int InnerCurrentControl<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT* y, const ScalarT* yp, const ScalarT* input, const ScalarT*, ScalarT* f)
      {
        const ScalarT vfd   = y[6];
        const ScalarT vfq   = y[7];
        const ScalarT icmdd = input[4] - input[6] * C_ * vfq;
        const ScalarT icmdq = input[5] + input[6] * C_ * vfd;
        const ScalarT li    = std::sqrt(Math::max(ONE<RealT>, ai_ * (icmdd * icmdd + icmdq * icmdq)));
        const ScalarT ed    = y[2] - input[2];
        const ScalarT eq    = y[3] - input[3];
        f[0]                = -yp[0] + Ki_ * ed + Kaw_ * (input[7] - y[4]);
        f[1]                = -yp[1] + Ki_ * eq + Kaw_ * (input[8] - y[5]);
        f[2]                = y[2] - icmdd / li;
        f[3]                = y[3] - icmdq / li;
        f[4]                = y[4] - (input[0] - input[6] * L_ * input[3] + Kp_ * ed + y[0]);
        f[5]                = y[5] - (input[1] + input[6] * L_ * input[2] + Kp_ * eq + y[1]);
        f[6]                = -yp[6] + (input[0] - vfd) / Tf_;
        f[7]                = -yp[7] + (input[1] - vfq) / Tf_;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int InnerCurrentControl<scalar_type, index_type>::evaluateInternalResidual()
      {
        this->gatherExternalVariables();
        evaluateInternalResidual(y_.getData(), yp_.getData(), y_ext_.data(), yp_ext_.data(), f_.getData());
        f_.setDataUpdated();
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int InnerCurrentControl<scalar_type, index_type>::evaluateResidual()
      {
        return evaluateInternalResidual();
      }

      template <typename scalar_type, typename index_type>
      void InnerCurrentControl<scalar_type, index_type>::initializeMonitor()
      {
        using Mon = typename ModelDataT::MonitorableVariables;
        monitor_->set(Mon::xid, [this]
                      { return y_.getData()[0]; });
        monitor_->set(Mon::xiq, [this]
                      { return y_.getData()[1]; });
        monitor_->set(Mon::ilimd, [this]
                      { return y_.getData()[2]; });
        monitor_->set(Mon::ilimq, [this]
                      { return y_.getData()[3]; });
        monitor_->set(Mon::ud, [this]
                      { return y_.getData()[4]; });
        monitor_->set(Mon::uq, [this]
                      { return y_.getData()[5]; });
      }

      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* InnerCurrentControl<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
