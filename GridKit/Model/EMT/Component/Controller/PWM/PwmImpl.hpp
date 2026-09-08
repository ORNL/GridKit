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
                      { return output(0); });
        monitor_->set(ModelDataT::MonitorableVariables::sb, [this]
                      { return output(1); });
        monitor_->set(ModelDataT::MonitorableVariables::sc, [this]
                      { return output(2); });
        for (size_t phase = 0; phase < 3; ++phase)
        {
          output_port_[phase].setComputed(
              [this, phase]
              { return output(phase); },
              [](typename SignalT::GradientT&, RealT) {});
        }
      }

      template <typename scalar_type, typename index_type>
      Pwm<scalar_type, index_type>::~Pwm() = default;

      template <typename scalar_type, typename index_type>
      void Pwm<scalar_type, index_type>::assignInput(size_t phase, SignalT* signal)
      {
        if (this->allocated_)
          throw std::logic_error("Assign PWM inputs before allocation");
        auto& input = input_.at(phase);
        if (signal == nullptr || (input != nullptr && input != signal))
          throw std::invalid_argument("Invalid PWM input assignment");
        input = signal;
      }

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
        this->allocated_ = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Pwm<scalar_type, index_type>::initialize(const std::map<Outputs, RealT>& outputs)
      {
        this->validateOutputValues(outputs);
        if (const int status = verify(); status != 0)
          return status;
        resetHistory();
        for (const auto& [key, value] : outputs)
        {
          this->checkOutputValue(outputs, key, static_cast<RealT>(output(static_cast<size_t>(key))));
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
        const RealT ratio            = fc_ / fm_;
        const RealT triple           = std::round(ratio / 3);
        const RealT tolerance        = 100 * std::numeric_limits<RealT>::epsilon() * ratio;
        sinusoidal_parameters_valid_ = parameters_valid_ && std::isfinite(M_) && M_ >= 0 && M_ <= 1
                                       && std::isfinite(fm_) && fm_ > 0 && fc_ > fm_
                                       && std::isfinite(ratio) && triple >= 1
                                       && std::abs(ratio - 3 * triple) <= tolerance
                                       && std::isfinite(1 / fm_)
                                       && ratio + horizon_ * fc_ < static_cast<RealT>(std::numeric_limits<long long>::max() / 2);
      }

      template <typename scalar_type, typename index_type>
      int Pwm<scalar_type, index_type>::verify() const
      {
        if (!parameters_valid_)
        {
          Log::error() << "PWM: require finite fc > 0 and 0 <= alignment <= 1\n";
          return 1;
        }
        if (sampledInput())
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
          Log::error() << "PWM: without modulation inputs require finite 0 <= M <= 1, fc > fm > 0, and fc/fm in 3N\n";
          return 1;
        }
        return 0;
      }

      template <typename scalar_type, typename index_type>
      bool Pwm<scalar_type, index_type>::sampledInput() const
      {
        return std::any_of(input_.begin(), input_.end(), [](const auto* signal)
                           { return signal != nullptr; });
      }

      template <typename scalar_type, typename index_type>
      auto Pwm<scalar_type, index_type>::readModulation() const -> std::array<RealT, 3>
      {
        std::array<RealT, 3> modulation;
        for (size_t phase = 0; phase < 3; ++phase)
        {
          const RealT value = static_cast<RealT>(input_[phase]->read());
          if (!std::isfinite(value) || std::abs(value) > 1 + 100 * std::numeric_limits<RealT>::epsilon())
            throw std::domain_error("PWM: sampled modulation must lie in [-1,1]");
          modulation[phase] = std::clamp(value, RealT{-1}, RealT{1});
        }
        return modulation;
      }

      template <typename scalar_type, typename index_type>
      void Pwm<scalar_type, index_type>::invalidateCache()
      {
        cached_time_.fill(std::numeric_limits<RealT>::quiet_NaN());
      }

      template <typename scalar_type, typename index_type>
      void Pwm<scalar_type, index_type>::resetHistory()
      {
        samples_.clear();
        invalidateCache();
        if (sampledInput())
        {
          if (verify() != 0 || !std::isfinite(this->time_)
              || std::abs(this->time_ * fc_) >= static_cast<RealT>(std::numeric_limits<long long>::max() / 2))
            throw std::domain_error("Cannot initialize sampled PWM with invalid inputs or time");
          // Initialize the active command and the one queued for the next interval.
          const auto interval   = static_cast<long long>(std::floor(this->time_ * fc_));
          const auto modulation = readModulation();
          samples_.push_back({interval, modulation});
          samples_.push_back({interval + 1, modulation});
        }
      }

      template <typename scalar_type, typename index_type>
      void Pwm<scalar_type, index_type>::acceptStep(RealT time)
      {
        if (!sampledInput())
          return;
        if (samples_.empty())
          throw std::logic_error("PWM: initialize held modulation before accepting steps");
        // Interval k+1 is sampled at t_k: one carrier of computational delay.
        const auto  next_interval = samples_.back().interval + 1;
        const RealT next          = static_cast<RealT>(next_interval - 1) / fc_;
        if (time < next)
          return;
        const RealT tolerance = 16 * std::numeric_limits<RealT>::epsilon() * std::max(RealT{1}, std::abs(next));
        if (time > next + tolerance)
          throw std::logic_error("PWM: solver skipped a modulation sampling instant");
        samples_.push_back({next_interval, readModulation()});
        // Retain only the active command and the one queued for the next interval.
        while (samples_.size() > 2)
          samples_.pop_front();
        invalidateCache();
      }

      template <typename scalar_type, typename index_type>
      auto Pwm<scalar_type, index_type>::nextDiscontinuityTime(RealT after) const -> RealT
      {
        if (!sampledInput())
          return std::numeric_limits<RealT>::infinity();
        if (samples_.empty())
          throw std::logic_error("PWM: initialize held modulation before requesting sampling times");
        const RealT next = static_cast<RealT>(samples_.back().interval) / fc_;
        if (next <= after)
          throw std::logic_error("PWM: uncommitted modulation sampling instant");
        return next;
      }

      template <typename scalar_type, typename index_type>
      auto Pwm<scalar_type, index_type>::maximumStepSize() const -> RealT
      {
        return sampledInput() ? std::min(1 / (20 * fc_), 1 / Math::MU<RealT>)
                              : std::numeric_limits<RealT>::infinity();
      }

      template <typename scalar_type, typename index_type>
      auto Pwm<scalar_type, index_type>::output(size_t phase) const -> ScalarT
      {
        if (phase >= 3 || verify() != 0 || !std::isfinite(this->time_))
        {
          throw std::domain_error("Cannot evaluate PWM with invalid parameters, time, or phase");
        }
        if (cached_time_[phase] == this->time_)
        {
          return static_cast<ScalarT>(cached_output_[phase]);
        }
        const RealT                pi = std::numbers::pi_v<RealT>;
        const std::array<RealT, 3> phi{0, -2 * pi / 3, 2 * pi / 3};
        const RealT                tc      = 1 / fc_;
        const bool                 sampled = sampledInput();
        if (sampled && samples_.empty())
          throw std::logic_error("PWM: initialize held modulation before evaluating switching outputs");
        const RealT t          = sampled ? this->time_ : std::remainder(this->time_, 1 / fm_);
        const auto  first      = static_cast<long long>(std::floor((t - horizon_) * fc_));
        const auto  last       = static_cast<long long>(std::floor((t + horizon_) * fc_));
        const auto  intervals  = sampled ? 0LL : static_cast<long long>(std::round(fc_ / fm_));
        RealT       sum        = 0;
        RealT       correction = 0;
        for (auto k = first; k <= last; ++k)
        {
          RealT modulation;
          if (sampled)
          {
            // Every periodic pulse replica uses the active held command.
            modulation = samples_.front().modulation[phase];
          }
          else
          {
            // Reduce the prescribed sinusoid while retaining carrier prehistory.
            const auto sample = static_cast<RealT>(k % intervals) + alignment_;
            modulation        = M_ * std::sin(2 * pi * sample / static_cast<RealT>(intervals) + phi[phase]);
          }
          const RealT duty       = (1 + modulation) / 2;
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
