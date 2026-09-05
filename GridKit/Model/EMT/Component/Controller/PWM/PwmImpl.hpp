#pragma once

#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

#include <GridKit/Model/EMT/Component/Controller/PWM/Pwm.hpp>
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
                      { return output(0); });
        monitor_->set(ModelDataT::MonitorableVariables::sb, [this]
                      { return output(1); });
        monitor_->set(ModelDataT::MonitorableVariables::sc, [this]
                      { return output(2); });
        for (size_t phase = 0; phase < 3; ++phase)
        {
          output_port_.signals[phase].setComputed(
              [this, phase]
              { return output(phase); },
              [](typename SignalT::GradientT&, RealT) {});
        }
      }

      template <typename scalar_type, typename index_type>
      Pwm<scalar_type, index_type>::~Pwm() = default;

      template <typename scalar_type, typename index_type>
      void Pwm<scalar_type, index_type>::assignOutput(size_t phase, SignalT* signal)
      {
        auto& assigned = assigned_output_.at(phase);
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
            { return output_port_.signals[phase].read(); },
            [this, phase](typename SignalT::GradientT& gradient, RealT scale)
            { output_port_.signals[phase].appendGradient(gradient, scale); });
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
        this->allocated_ = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Pwm<scalar_type, index_type>::initialize()
      {
        return verify();
      }

      template <typename scalar_type, typename index_type>
      int Pwm<scalar_type, index_type>::tagDifferentiable()
      {
        return 0;
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
      int Pwm<scalar_type, index_type>::evaluateJacobian()
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
        using Parameter = typename ModelDataT::Parameters;
        auto read       = [&data](Parameter key, RealT fallback)
        {
          const auto entry = data.parameters.find(key);
          if (entry == data.parameters.end())
          {
            return fallback;
          }
          if (const auto* value = std::get_if<RealT>(&entry->second))
          {
            return *value;
          }
          if (const auto* value = std::get_if<IdxT>(&entry->second))
          {
            return static_cast<RealT>(*value);
          }
          return std::numeric_limits<RealT>::quiet_NaN();
        };
        const auto missing = std::numeric_limits<RealT>::quiet_NaN();
        M_                 = read(Parameter::M, missing);
        fm_                = read(Parameter::fm, missing);
        fc_                = read(Parameter::fc, missing);
        alignment_         = read(Parameter::alignment, RealT{0.5});
        horizon_           = std::log(4 / std::numeric_limits<RealT>::epsilon()) / Math::MU<RealT>;

        const RealT ratio     = fc_ / fm_;
        const RealT triple    = std::round(ratio / 3);
        const RealT tolerance = 100 * std::numeric_limits<RealT>::epsilon() * ratio;
        parameters_valid_     = std::isfinite(M_) && M_ >= 0 && M_ <= 1
                            && std::isfinite(fm_) && fm_ > 0
                            && std::isfinite(fc_) && fc_ > fm_
                            && std::isfinite(alignment_) && alignment_ >= 0 && alignment_ <= 1
                            && std::isfinite(ratio) && triple >= 1
                            && std::abs(ratio - 3 * triple) <= tolerance
                            && std::isfinite(1 / fm_) && std::isfinite(1 / fc_)
                            && ratio + horizon_ * fc_ < static_cast<RealT>(std::numeric_limits<long long>::max() / 2);
      }

      template <typename scalar_type, typename index_type>
      int Pwm<scalar_type, index_type>::verify() const
      {
        if (!parameters_valid_)
        {
          Log::error() << "PWM: require finite 0 <= M <= 1, fc > fm > 0, fc/fm in 3N, and 0 <= alignment <= 1\n";
          return 1;
        }
        return 0;
      }

      template <typename scalar_type, typename index_type>
      auto Pwm<scalar_type, index_type>::output(size_t phase) const -> ScalarT
      {
        if (phase >= 3 || !parameters_valid_ || !std::isfinite(this->time_))
        {
          throw std::domain_error("Cannot evaluate PWM with invalid parameters, time, or phase");
        }
        if (cached_time_[phase] == this->time_)
        {
          return static_cast<ScalarT>(cached_output_[phase]);
        }
        const RealT                pi = std::numbers::pi_v<RealT>;
        const std::array<RealT, 3> phi{0, -2 * pi / 3, 2 * pi / 3};
        const RealT                tc         = 1 / fc_;
        const RealT                t          = std::remainder(this->time_, 1 / fm_);
        const auto                 first      = static_cast<long long>(std::floor((t - horizon_) * fc_));
        const auto                 last       = static_cast<long long>(std::floor((t + horizon_) * fc_));
        const auto                 intervals  = static_cast<long long>(std::round(fc_ / fm_));
        RealT                      sum        = 0;
        RealT                      correction = 0;
        for (auto k = first; k <= last; ++k)
        {
          // Reduce the modulation argument while retaining the carrier interval
          // on the whole time axis, including the prehistory before t = 0.
          const auto  sample     = static_cast<RealT>(k % intervals) + alignment_;
          const RealT duty       = (1 + M_ * std::sin(2 * pi * sample / static_cast<RealT>(intervals) + phi[phase])) / 2;
          const RealT local_time = t - static_cast<RealT>(k) * tc;
          const RealT on         = alignment_ * (1 - duty) * tc;
          const RealT off        = (alignment_ + (1 - alignment_) * duty) * tc;
          // Reflection avoids subtracting two values close to one in the tail.
          const RealT pulse      = local_time <= (on + off) / 2
                                       ? Math::sigmoid(local_time - on) - Math::sigmoid(local_time - off)
                                       : Math::sigmoid(off - local_time) - Math::sigmoid(on - local_time);
          const RealT term       = pulse - correction;
          const RealT next       = sum + term;
          correction             = (next - sum) - term;
          sum                    = next;
        }
        // The omitted tails are bounded by 2 sigmoid(-horizon), independently
        // of the carrier frequency, since the unsmoothed pulses do not overlap.
        cached_time_[phase]   = this->time_;
        cached_output_[phase] = std::clamp(sum, RealT{0}, RealT{1});
        return static_cast<ScalarT>(cached_output_[phase]);
      }
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
