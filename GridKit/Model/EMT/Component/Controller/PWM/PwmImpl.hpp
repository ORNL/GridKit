#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Component/Controller/PWM/Pwm.hpp>
#include <GridKit/Model/EMT/ComponentInitialization.hpp>
#include <GridKit/Model/EMT/Operators/Reference/Park/ParkImpl.hpp>
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
        using Mon = typename ModelDataT::MonitorableVariables;
        monitor_->set(Mon::sa, [this]
                      { return output_port_[0].read(); });
        monitor_->set(Mon::sb, [this]
                      { return output_port_[1].read(); });
        monitor_->set(Mon::sc, [this]
                      { return output_port_[2].read(); });
        monitor_->set(Mon::ma, [this]
                      { return modulation(0); });
        monitor_->set(Mon::mb, [this]
                      { return modulation(1); });
        monitor_->set(Mon::mc, [this]
                      { return modulation(2); });
        monitor_->set(Mon::ulimd, [this]
                      { return output_port_[3].read(); });
        monitor_->set(Mon::ulimq, [this]
                      { return output_port_[4].read(); });
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
      void Pwm<scalar_type, index_type>::attachInput(
          const std::array<SignalT*, 2>& command, SignalT* vdc, SignalT* theta)
      {
        if (this->allocated_)
          throw std::logic_error("Attach PWM inputs before allocation");
        input_ = {command[0], command[1], vdc, theta};
      }

      template <typename scalar_type, typename index_type>
      void Pwm<scalar_type, index_type>::assignOutput(Outputs output, SignalT* signal)
      {
        const auto index    = static_cast<size_t>(output);
        auto&      assigned = assigned_output_.at(index);
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
            [this, index]
            { return output_port_[index].read(); },
            [this, index](typename SignalT::GradientT& gradient, RealT scale)
            { output_port_[index].appendGradient(gradient, scale); });
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
        Mmax_              = parameter<RealT>(data, Parameter::Mmax, RealT{1});
        if (Mmax_ > ZERO<RealT>)
        {
          au_ = RealT{8} / (RealT{3} * Mmax_ * Mmax_);
        }

        parameters_valid_ = std::isfinite(fc_) && fc_ > 0
                            && std::isfinite(1 / fc_)
                            && std::isfinite(alignment_) && alignment_ >= 0 && alignment_ <= 1
                            && std::isfinite(Mmax_) && Mmax_ > 0 && Mmax_ <= 1
                            && std::isfinite(au_) && au_ > 0
                            && std::isfinite(Math::MU<RealT>) && Math::MU<RealT> > 0;
        sinusoidal_parameters_valid_ = parameters_valid_ && std::isfinite(M_) && M_ >= 0 && M_ <= 1
                                       && std::isfinite(fm_) && fm_ > 0 && fc_ > fm_
                                       && std::isfinite(1 / fm_);
      }

      template <typename scalar_type, typename index_type>
      int Pwm<scalar_type, index_type>::verify() const
      {
        if (!parameters_valid_)
        {
          Log::error() << "PWM: require finite fc > 0, 0 <= alignment <= 1, and 0 < Mmax <= 1\n";
          return 1;
        }
        if (hasInput())
        {
          for (const auto* signal : input_)
            if (signal == nullptr || !signal->linked())
            {
              Log::error() << "PWM: voltage command, DC voltage, and angle inputs must all have linked sources\n";
              return 1;
            }
          return 0;
        }
        if (!sinusoidal_parameters_valid_)
        {
          Log::error() << "PWM: without voltage command inputs require finite 0 <= M <= 1 and fc > fm > 0\n";
          return 1;
        }
        if (assigned_output_[3] != nullptr || assigned_output_[4] != nullptr)
        {
          Log::error() << "PWM: the limited voltage command requires the voltage command inputs\n";
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
      auto Pwm<scalar_type, index_type>::inputValues() const -> std::array<ScalarT, 4>
      {
        std::array<ScalarT, 4> values{};
        for (size_t n = 0; n < values.size(); ++n)
        {
          values[n] = input_[n]->read();
          if (!std::isfinite(static_cast<RealT>(values[n])))
            throw std::domain_error("PWM: voltage command, DC voltage, and angle must be finite");
        }
        if (static_cast<RealT>(values[2]) < ZERO<RealT>)
          throw std::domain_error("PWM: DC voltage must be nonnegative");
        return values;
      }

      template <typename scalar_type, typename index_type>
      auto Pwm<scalar_type, index_type>::fraction(const ScalarT* input) const -> std::array<ScalarT, 2>
      {
        const ScalarT scale = std::sqrt(Math::max(input[2] * input[2], au_ * (input[0] * input[0] + input[1] * input[1])));
        return {input[0] / scale, input[1] / scale};
      }

      template <typename scalar_type, typename index_type>
      auto Pwm<scalar_type, index_type>::phaseModulation(size_t phase, const ScalarT* input) const -> ScalarT
      {
        const auto w      = fraction(input);
        const auto matrix = Park<ScalarT, IdxT>::template transformation<ScalarT>(input[3]);
        return ScalarT{2} * (matrix[0][phase] * w[0] + matrix[1][phase] * w[1]);
      }

      template <typename scalar_type, typename index_type>
      auto Pwm<scalar_type, index_type>::modulation(size_t phase) const -> ScalarT
      {
        if (phase >= 3 || verify() != 0 || !std::isfinite(this->time_))
          throw std::domain_error("PWM: invalid modulation phase, configuration, or time");
        if (hasInput())
          return phaseModulation(phase, inputValues().data());
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
      auto Pwm<scalar_type, index_type>::switching(ScalarT duty) const -> ScalarT
      {
        const RealT resolution = Math::MU<RealT> / fc_;
        // The omitted periodic ripple is below double roundoff.
        if (resolution < RealT{0.4})
          return duty;
        const int   radius = resolution >= RealT{16} ? 4 : (resolution >= RealT{4} ? 16 : 127);
        const RealT tc     = 1 / fc_;
        const RealT t      = std::remainder(this->time_, tc);
        ScalarT     sum{0};
        for (int k = -radius; k <= radius; ++k)
          sum += pulse(duty, t - k * tc);
        return sum;
      }

      template <typename scalar_type, typename index_type>
      auto Pwm<scalar_type, index_type>::evaluateOutput(Outputs output, const ScalarT* input) const -> ScalarT
      {
        const auto index = static_cast<size_t>(output);
        if (index >= 3)
          return input[2] * fraction(input)[index - 3];
        return switching((ScalarT{1} + phaseModulation(index, input)) / 2);
      }

      template <typename scalar_type, typename index_type>
      auto Pwm<scalar_type, index_type>::output(Outputs output) const -> ScalarT
      {
        const auto index = static_cast<size_t>(output);
        if (index >= output_port_.size() || verify() != 0 || !std::isfinite(this->time_))
          throw std::domain_error("PWM: invalid output, configuration, or time");
        if (hasInput())
          return evaluateOutput(output, inputValues().data());
        if (index >= 3)
          throw std::logic_error("PWM: limited voltage requires voltage command inputs");
        return switching((ScalarT{1} + modulation(index)) / 2);
      }
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
