#pragma once

#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

#include <GridKit/Model/EMT/Component/Controller/PWM/Pwm.hpp>
#include <GridKit/Model/EMT/ComponentInitialization.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template <typename scalar_type, typename index_type>
      Pwm<scalar_type, index_type>::Pwm()
        : Pwm(ModelDataT{})
      {
      }

      template <typename scalar_type, typename index_type>
      Pwm<scalar_type, index_type>::Pwm(const ModelDataT& data)
        : monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        monitor_->set(ModelDataT::MonitorableVariables::sa, [this]
                      { return output(Outputs::sa); });
        monitor_->set(ModelDataT::MonitorableVariables::sb, [this]
                      { return output(Outputs::sb); });
        monitor_->set(ModelDataT::MonitorableVariables::sc, [this]
                      { return output(Outputs::sc); });
        for (size_t n = 0; n < output_port_.size(); ++n)
        {
          const auto key = static_cast<Outputs>(n);
          output_port_[n].setComputed(
              [this, key]
              { return output(key); },
              [this, key](typename SignalT::GradientT& gradient, RealT scale)
              { appendOutputGradient(key, gradient, scale); });
        }
      }

      template <typename scalar_type, typename index_type>
      Pwm<scalar_type, index_type>::~Pwm() = default;

      template <typename scalar_type, typename index_type>
      void Pwm<scalar_type, index_type>::assignInput(Inputs key, SignalT* signal)
      {
        if (this->allocated_)
          throw std::logic_error("Assign PWM inputs before allocation");
        auto& input = input_.at(static_cast<size_t>(key));
        if (signal == nullptr || (input != nullptr && input != signal))
          throw std::invalid_argument("Invalid PWM input assignment");
        input = signal;
      }

      template <typename scalar_type, typename index_type>
      void Pwm<scalar_type, index_type>::assignOutput(Outputs output, SignalT* signal)
      {
        const auto phase    = static_cast<size_t>(output);
        auto&      assigned = assigned_output_.at(phase);
        if (signal == nullptr || (assigned != nullptr && assigned != signal))
        {
          throw std::invalid_argument("Invalid Pwm output assignment");
        }
        if (assigned == signal)
        {
          return;
        }
        signal->claimProducer();
        assigned = signal;
        signal->setComputed(
            [this, phase]
            { return output_port_[phase].read(); },
            [this, phase](typename SignalT::GradientT& gradient, RealT scale)
            { output_port_[phase].appendGradient(gradient, scale); });
      }

      template <typename scalar_type, typename index_type>
      int Pwm<scalar_type, index_type>::setGridKitComponentID(IdxT id)
      {
        this->gridkit_component_id_ = id;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Pwm<scalar_type, index_type>::allocate()
      {
        if (hasInput())
        {
          this->allocateExternalVectors(static_cast<IdxT>(input_.size()), 0);
          for (IdxT n = 0; n < static_cast<IdxT>(input_.size()); ++n)
          {
            this->setExternalVariableSignal(n, input_[static_cast<size_t>(n)]);
          }
        }
        this->allocated_ = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Pwm<scalar_type, index_type>::initialize(const std::map<Outputs, RealT>& outputs)
      {
        this->validateOutputValues(outputs);
        if (const int status = verify(); status != 0)
          return status;
        for (const auto& [key, value] : outputs)
        {
          this->checkOutputValue(outputs, key, static_cast<RealT>(output(key)));
        }
        return verify();
      }

      template <typename scalar_type, typename index_type>
      int Pwm<scalar_type, index_type>::setAbsoluteTolerance(RealT)
      {
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Pwm<scalar_type, index_type>::evaluateInternalResidual()
      {
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Pwm<scalar_type, index_type>::evaluateExternalResidual()
      {
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Pwm<scalar_type, index_type>::evaluateResidual()
      {
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Pwm<scalar_type, index_type>::assembleJacobian(RealT, RealT)
      {
        return 0;
      }

      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* Pwm<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      template <typename scalar_type, typename index_type>
      void Pwm<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Parameter    = typename ModelDataT::Parameters;
        const auto missing = std::numeric_limits<RealT>::quiet_NaN();
        M_                 = parameter<RealT>(data, Parameter::M, missing);
        fm_                = parameter<RealT>(data, Parameter::fm, missing);
        fc_                = parameter<RealT>(data, Parameter::fc, missing);
        alignment_         = parameter<RealT>(data, Parameter::alignment, RealT{0.5});
        horizon_           = std::log(4 / std::numeric_limits<RealT>::epsilon()) / Math::MU<RealT>;

        parameters_valid_ = std::isfinite(fc_) && fc_ > 0
                            && std::isfinite(1 / fc_)
                            && std::isfinite(alignment_) && alignment_ >= 0 && alignment_ <= 1
                            && std::isfinite(horizon_) && horizon_ > 0
                            && horizon_ * fc_ < static_cast<RealT>(std::numeric_limits<long long>::max() / 2);
        sinusoidal_parameters_valid_ = parameters_valid_ && std::isfinite(M_) && M_ >= 0 && M_ <= 1
                                       && std::isfinite(fm_) && fm_ > 0 && fc_ > fm_
                                       && std::isfinite(1 / fm_);
      }

      template <typename scalar_type, typename index_type>
      int Pwm<scalar_type, index_type>::verify() const
      {
        if (!parameters_valid_)
        {
          Log::error() << "PWM: require finite fc > 0 and 0 <= alignment <= 1\n";
          return 1;
        }
        if (hasInput())
        {
          for (const auto* signal : input_)
            if (signal == nullptr || !signal->linked())
            {
              Log::error() << "PWM: all three modulation inputs must have linked sources\n";
              return 1;
            }
        }
        else if (!sinusoidal_parameters_valid_)
        {
          Log::error() << "PWM: without modulation inputs require finite 0 <= M <= 1 and fc > fm > 0\n";
          return 1;
        }
        return 0;
      }

      template <typename scalar_type, typename index_type>
      bool Pwm<scalar_type, index_type>::hasInput() const
      {
        return std::any_of(input_.begin(), input_.end(), [](const auto* signal)
                           { return signal != nullptr; });
      }

      template <typename scalar_type, typename index_type>
      auto Pwm<scalar_type, index_type>::modulation(size_t phase) const -> ScalarT
      {
        if (hasInput())
        {
          const auto value = input_[phase]->read();
          if (!std::isfinite(static_cast<RealT>(value)) || std::abs(static_cast<RealT>(value)) > 1)
            throw std::domain_error("PWM: modulation must lie in [-1,1]");
          return value;
        }
        const RealT                pi = std::numbers::pi_v<RealT>;
        const std::array<RealT, 3> phi{0, -2 * pi / 3, 2 * pi / 3};
        return ScalarT{M_ * std::sin(2 * pi * fm_ * std::remainder(this->time_, 1 / fm_) + phi[phase])};
      }

      template <typename scalar_type, typename index_type>
      auto Pwm<scalar_type, index_type>::maximumStepSize() const -> RealT
      {
        return 1 / Math::MU<RealT>;
      }

      template <typename scalar_type, typename index_type>
      auto Pwm<scalar_type, index_type>::pulse(ScalarT duty, RealT local_time) const -> ScalarT
      {
        const RealT tc  = 1 / fc_;
        const auto  on  = alignment_ * (1 - duty) * tc;
        const auto  off = (alignment_ + (1 - alignment_) * duty) * tc;
        // Reflection avoids subtracting two values close to one in the tail.
        if (local_time <= static_cast<RealT>((on + off) / 2))
          return Math::sigmoid(local_time - on) - Math::sigmoid(local_time - off);
        return Math::sigmoid(off - local_time) - Math::sigmoid(on - local_time);
      }

      template <typename scalar_type, typename index_type>
      auto Pwm<scalar_type, index_type>::output(Outputs output) const -> ScalarT
      {
        const auto phase = static_cast<size_t>(output);
        if (phase >= 3 || verify() != 0 || !std::isfinite(this->time_))
        {
          throw std::domain_error("Cannot evaluate PWM with invalid parameters, time, or phase");
        }
        const auto  duty  = (1 + modulation(phase)) / 2;
        const RealT tc    = 1 / fc_;
        const RealT t     = std::remainder(this->time_, tc);
        const auto  first = static_cast<long long>(std::floor((t - horizon_) * fc_));
        const auto  last  = static_cast<long long>(std::floor((t + horizon_) * fc_));
        ScalarT     sum{0};
        ScalarT     correction{0};
        for (auto k = first; k <= last; ++k)
        {
          const auto term = pulse(duty, t - static_cast<RealT>(k) * tc) - correction;
          const auto next = sum + term;
          correction      = (next - sum) - term;
          sum             = next;
        }
        return sum;
      }

      template <typename scalar_type, typename index_type>
      void Pwm<scalar_type, index_type>::appendOutputGradient(
          Outputs output, typename SignalT::GradientT& gradient, RealT scale) const
      {
        const auto phase = static_cast<size_t>(output);
        if (phase >= 3 || verify() != 0 || !std::isfinite(this->time_))
          throw std::domain_error("Cannot differentiate PWM with invalid parameters, time, or phase");
        if (!hasInput())
          return;
        const RealT duty  = (1 + static_cast<RealT>(modulation(phase))) / 2;
        const RealT tc    = 1 / fc_;
        const RealT t     = std::remainder(this->time_, tc);
        const RealT on    = alignment_ * (1 - duty) * tc;
        const RealT off   = (alignment_ + (1 - alignment_) * duty) * tc;
        const auto  first = static_cast<long long>(std::floor((t - horizon_) * fc_));
        const auto  last  = static_cast<long long>(std::floor((t + horizon_) * fc_));
        RealT       derivative{0};
        for (auto k = first; k <= last; ++k)
        {
          const RealT local  = t - static_cast<RealT>(k) * tc;
          const RealT a      = Math::sigmoid(-std::abs(local - on));
          const RealT b      = Math::sigmoid(-std::abs(local - off));
          derivative        += alignment_ * a * (1 - a) + (1 - alignment_) * b * (1 - b);
        }
        input_[phase]->appendGradient(gradient, scale * tc * Math::MU<RealT> * derivative / 2);
      }
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
