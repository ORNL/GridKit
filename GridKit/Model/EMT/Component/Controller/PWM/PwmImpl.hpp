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
        horizon_           = std::log(4 / std::numeric_limits<RealT>::epsilon()) / Math::MU<RealT>;
        if (Mmax_ > ZERO<RealT>)
        {
          au_ = RealT{8} / (RealT{3} * Mmax_ * Mmax_);
        }

        parameters_valid_ = std::isfinite(fc_) && fc_ > 0
                            && std::isfinite(1 / fc_)
                            && std::isfinite(alignment_) && alignment_ >= 0 && alignment_ <= 1
                            && std::isfinite(Mmax_) && Mmax_ > 0 && Mmax_ <= 1
                            && std::isfinite(horizon_) && horizon_ > 0
                            && horizon_ * fc_ < static_cast<RealT>(std::numeric_limits<long long>::max() / 2);
        sinusoidal_parameters_valid_ = parameters_valid_ && std::isfinite(M_) && M_ >= 0 && M_ <= 1
                                       && std::isfinite(fm_) && fm_ > 0 && fc_ > fm_
                                       && std::isfinite(1 / fm_);
        if (parameters_valid_)
          replica_decay_ = std::exp(-Math::MU<RealT> / fc_);
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
      auto Pwm<scalar_type, index_type>::fraction() const -> std::array<ScalarT, 2>
      {
        const ScalarT vdc = input_[2]->read();
        if (!std::isfinite(static_cast<RealT>(vdc)) || static_cast<RealT>(vdc) < ZERO<RealT>)
          throw std::domain_error("PWM: DC voltage must be finite and nonnegative");
        const ScalarT ud = input_[0]->read();
        const ScalarT uq = input_[1]->read();
        if (!std::isfinite(static_cast<RealT>(ud)) || !std::isfinite(static_cast<RealT>(uq))
            || !std::isfinite(static_cast<RealT>(input_[3]->read())))
          throw std::domain_error("PWM: voltage command and angle inputs must be finite");
        // Radial limit of the command against the available DC voltage.
        const ScalarT scale = std::sqrt(Math::max(vdc * vdc, au_ * (ud * ud + uq * uq)));
        return {ud / scale, uq / scale};
      }

      template <typename scalar_type, typename index_type>
      auto Pwm<scalar_type, index_type>::modulation(size_t phase) const -> ScalarT
      {
        if (hasInput())
        {
          const auto w      = fraction();
          const auto matrix = Park<ScalarT, IdxT>::template transformation<ScalarT>(input_[3]->read());
          return ScalarT{2} * (matrix[0][phase] * w[0] + matrix[1][phase] * w[1]);
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
        const auto index = static_cast<size_t>(output);
        if (index >= output_port_.size() || verify() != 0 || !std::isfinite(this->time_))
        {
          throw std::domain_error("Cannot evaluate PWM with invalid parameters, time, or output");
        }
        if (index >= 3)
        {
          if (!hasInput())
            throw std::domain_error("PWM: the limited voltage command requires the voltage command inputs");
          return input_[2]->read() * fraction()[index - 3];
        }
        const auto  duty  = (1 + modulation(index)) / 2;
        const RealT tc    = 1 / fc_;
        const RealT t     = std::remainder(this->time_, tc);
        const auto  first = static_cast<long long>(std::floor((t - horizon_) * fc_));
        const auto  last  = static_cast<long long>(std::floor((t + horizon_) * fc_));
        const RealT mu    = Math::MU<RealT>;
        const auto  on    = alignment_ * (1 - duty) * tc;
        const auto  off   = (alignment_ + (1 - alignment_) * duty) * tc;

        // Evaluate the nearest pulse directly; the remaining replicas form two tails.
        const auto center = std::clamp(
            static_cast<long long>(std::round((t - static_cast<RealT>((on + off) / 2)) * fc_)),
            first,
            last);
        const auto width = mu * duty * tc;
        const auto r     = std::exp(-width);
        const auto h     = std::tanh(width / 2);
        const auto span  = 2 * h / (1 + h); // Stable 1 - exp(-width).

        ScalarT                        sum = pulse(duty, t - static_cast<RealT>(center) * tc);
        ScalarT                        correction{0};
        const std::array<long long, 2> count{center - first, last - center};
        const std::array<ScalarT, 2>   distance{
            t - static_cast<RealT>(center - 1) * tc - off,
            static_cast<RealT>(center + 1) * tc + on - t};

        for (size_t side = 0; side < count.size(); ++side)
        {
          if (count[side] == 0)
            continue;
          auto z = std::exp(-mu * distance[side]);
          for (long long k = 0; k < count[side]; ++k)
          {
            const auto term  = z * span / ((1 + z) * (1 + z * r)) - correction;
            const auto next  = sum + term;
            correction       = (next - sum) - term;
            sum              = next;
            z               *= replica_decay_;
          }
        }
        return sum;
      }

      template <typename scalar_type, typename index_type>
      void Pwm<scalar_type, index_type>::appendOutputGradient(
          Outputs output, typename SignalT::GradientT& gradient, RealT scale) const
      {
        const auto index = static_cast<size_t>(output);
        if (index >= output_port_.size() || verify() != 0 || !std::isfinite(this->time_))
          throw std::domain_error("Cannot differentiate PWM with invalid parameters, time, or output");
        if (!hasInput())
          return;
        // Limiter geometry: w = u / sqrt(max(vdc^2, a_u |u|^2)) with a logistic gate on the smooth maximum.
        const RealT                vdc = static_cast<RealT>(input_[2]->read());
        const std::array<RealT, 2> u{static_cast<RealT>(input_[0]->read()), static_cast<RealT>(input_[1]->read())};
        const RealT                radius = au_ * (u[0] * u[0] + u[1] * u[1]);
        const RealT                square = Math::max(vdc * vdc, radius);
        const RealT                norm   = std::sqrt(square);
        const RealT                gate   = Math::sigmoid(radius - vdc * vdc);
        const std::array<RealT, 2> w{u[0] / norm, u[1] / norm};
        const auto                 fraction_command = [&](size_t k, size_t j)
        {
          RealT partial = -w[k] * gate * au_ * u[j] / square;
          if (k == j)
          {
            partial += 1 / norm;
          }
          return partial;
        };
        const auto fraction_dc = [&](size_t k)
        {
          return -w[k] * (1 - gate) * vdc / square;
        };
        if (index >= 3)
        {
          const size_t k = index - 3;
          for (size_t j = 0; j < 2; ++j)
            input_[j]->appendGradient(gradient, scale * vdc * fraction_command(k, j));
          input_[2]->appendGradient(gradient, scale * (w[k] + vdc * fraction_dc(k)));
          return;
        }
        // Switching function: pulse slope, then the inverse Park transform and the limiter.
        const size_t phase = index;
        const RealT  duty  = (1 + static_cast<RealT>(modulation(phase))) / 2;
        const RealT  tc    = 1 / fc_;
        const RealT  t     = std::remainder(this->time_, tc);
        const RealT  on    = alignment_ * (1 - duty) * tc;
        const RealT  off   = (alignment_ + (1 - alignment_) * duty) * tc;
        const auto   first = static_cast<long long>(std::floor((t - horizon_) * fc_));
        const auto   last  = static_cast<long long>(std::floor((t + horizon_) * fc_));
        RealT        derivative{0};
        for (auto k = first; k <= last; ++k)
        {
          const RealT local  = t - static_cast<RealT>(k) * tc;
          const RealT a      = Math::sigmoid(-std::abs(local - on));
          const RealT b      = Math::sigmoid(-std::abs(local - off));
          derivative        += alignment_ * a * (1 - a) + (1 - alignment_) * b * (1 - b);
        }
        const RealT slope  = scale * tc * Math::MU<RealT> * derivative / 2;
        const auto  matrix = Park<ScalarT, IdxT>::template transformation<RealT>(static_cast<RealT>(input_[3]->read()));
        for (size_t j = 0; j < 2; ++j)
        {
          RealT partial{0};
          for (size_t k = 0; k < 2; ++k)
            partial += matrix[k][phase] * 2 * fraction_command(k, j);
          input_[j]->appendGradient(gradient, slope * partial);
        }
        RealT dc_partial{0};
        for (size_t k = 0; k < 2; ++k)
          dc_partial += matrix[k][phase] * 2 * fraction_dc(k);
        input_[2]->appendGradient(gradient, slope * dc_partial);
        // The angle derivative of the cosine row is the sine row and vice versa.
        input_[3]->appendGradient(gradient, slope * 2 * (matrix[1][phase] * w[0] - matrix[0][phase] * w[1]));
      }
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
