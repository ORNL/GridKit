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
        size_ = 6;
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
        L_              = parameter<RealT>(data, Parameter::L, L_);
        Kp_             = parameter<RealT>(data, Parameter::Kp, Kp_);
        Ki_             = parameter<RealT>(data, Parameter::Ki, Ki_);
        Kaw_            = parameter<RealT>(data, Parameter::Kaw, Kaw_);
        Imax_           = parameter<RealT>(data, Parameter::Imax, Imax_);
        Mmax_           = parameter<RealT>(data, Parameter::Mmax, Mmax_);
        if (Imax_ > ZERO<RealT> && Mmax_ > ZERO<RealT>)
        {
          ai_ = ONE<RealT> / (Imax_ * Imax_);
          au_ = RealT{8} / (RealT{3} * Mmax_ * Mmax_);
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
        const auto inputs = signals_.attachedSignals({V::VD, V::VQ, V::ID, V::IQ, V::ICMDD, V::ICMDQ, V::OMEGA, V::VDC});
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
        for (const auto value : {L_, Kp_, Ki_, Kaw_, Imax_, Mmax_, ai_, au_})
          if (!std::isfinite(value) || value <= ZERO<RealT>)
          {
            Log::error() << "InnerCurrentControl: parameters must be finite and positive\n";
            ++error_count;
          }
        if (Mmax_ > ONE<RealT>)
        {
          Log::error() << "InnerCurrentControl: modulation limit must not exceed one\n";
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
        const RealT icmdd = static_cast<RealT>(signals_.template readExternalVariable<V::ICMDD>());
        const RealT icmdq = static_cast<RealT>(signals_.template readExternalVariable<V::ICMDQ>());
        const RealT omega = static_cast<RealT>(signals_.template readExternalVariable<V::OMEGA>());
        const RealT vdc   = static_cast<RealT>(signals_.template readExternalVariable<V::VDC>());
        for (const auto value : {vd, vq, id, iq, icmdd, icmdq, omega, vdc})
          if (!std::isfinite(value))
            throw std::invalid_argument("InnerCurrentControl: nonfinite initial input");
        if (vdc < ZERO<RealT>)
          throw std::invalid_argument("InnerCurrentControl: initial DC voltage must be nonnegative");
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
        const RealT                lu = std::sqrt(Math::max(vdc * vdc, au_ * (base[0] * base[0] + base[1] * base[1])));
        std::array<RealT, 2>       u{vdc * base[0] / lu, vdc * base[1] / lu};
        std::array<RealT, 2>       z = base;
        if (outputs.contains(Outputs::ud) || outputs.contains(Outputs::uq))
        {
          u[0]             = this->outputValue(outputs, Outputs::ud, u[0]);
          u[1]             = this->outputValue(outputs, Outputs::uq, u[1]);
          const RealT norm = std::hypot(u[0], u[1]);
          if (vdc == ZERO<RealT>)
          {
            if (norm != ZERO<RealT>)
              throw std::invalid_argument("InnerCurrentControl: initial voltage command must be zero at zero DC voltage");
          }
          else
          {
            if (norm >= vdc / std::sqrt(au_))
              throw std::invalid_argument("InnerCurrentControl: initial voltage command must be inside the modulation limit");
            // Invert the radial smooth limiter to recover the integral contribution.
            // A clipped output alone cannot determine a unique integral state.
            const auto gain = [&](RealT scale)
            {
              const RealT radius = scale * norm;
              return vdc * scale / std::sqrt(Math::max(vdc * vdc, au_ * radius * radius));
            };
            RealT lower = ZERO<RealT>;
            RealT upper = ONE<RealT>;
            while (gain(upper) < ONE<RealT> - tolerance)
            {
              upper *= RealT{2};
              if (!std::isfinite(upper) || !std::isfinite(gain(upper)))
                throw std::invalid_argument("InnerCurrentControl: initial voltage command cannot be inverted");
            }
            for (size_t n = 0; n < 64; ++n)
            {
              const RealT middle = lower + (upper - lower) / RealT{2};
              if (gain(middle) < ONE<RealT>)
                lower = middle;
              else
                upper = middle;
            }
            const RealT scale = lower + (upper - lower) / RealT{2};
            z                 = {scale * u[0], scale * u[1]};
          }
        }
        const std::array<RealT, 6> values{z[0] - base[0], z[1] - base[1], limited[0], limited[1], u[0], u[1]};
        for (const RealT value : values)
          if (!std::isfinite(value))
            throw std::invalid_argument("InnerCurrentControl: nonfinite derived initial value");
        auto* y = y_.getData();
        for (size_t n = 0; n < values.size(); ++n)
          y[n] = values[n];
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
        return {signals_.attachedSignals({V::VD, V::VQ, V::ID, V::IQ, V::ICMDD, V::ICMDQ, V::OMEGA, V::VDC}), {}, {}};
      }

      template <typename scalar_type, typename index_type>
      int InnerCurrentControl<scalar_type, index_type>::setAbsoluteTolerance(RealT tolerance)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(tolerance));
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int InnerCurrentControl<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT* y, const ScalarT* yp, const ScalarT* input, const ScalarT*, ScalarT* f)
      {
        const ScalarT li = std::sqrt(Math::max(ONE<RealT>, ai_ * (input[4] * input[4] + input[5] * input[5])));
        const ScalarT ed = y[2] - input[2];
        const ScalarT eq = y[3] - input[3];
        const ScalarT zd = input[0] - input[6] * L_ * input[3] + Kp_ * ed + y[0];
        const ScalarT zq = input[1] + input[6] * L_ * input[2] + Kp_ * eq + y[1];
        const ScalarT lu = std::sqrt(Math::max(input[7] * input[7], au_ * (zd * zd + zq * zq)));
        f[0]             = -yp[0] + Ki_ * ed + Kaw_ * (y[4] - zd);
        f[1]             = -yp[1] + Ki_ * eq + Kaw_ * (y[5] - zq);
        f[2]             = y[2] - input[4] / li;
        f[3]             = y[3] - input[5] / li;
        f[4]             = y[4] - input[7] * zd / lu;
        f[5]             = y[5] - input[7] * zq / lu;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int InnerCurrentControl<scalar_type, index_type>::evaluateInternalResidual()
      {
        this->gatherExternalVariables();
        const RealT vdc = static_cast<RealT>(y_ext_[7]);
        if (!std::isfinite(vdc) || vdc < ZERO<RealT>)
          throw std::invalid_argument("InnerCurrentControl: DC voltage must be finite and nonnegative");
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
