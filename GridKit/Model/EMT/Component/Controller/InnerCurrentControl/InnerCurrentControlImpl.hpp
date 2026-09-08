#pragma once

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Component/Controller/InnerCurrentControl/InnerCurrentControl.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template <typename scalar_type, typename index_type>
      InnerCurrentControl<scalar_type, index_type>::InnerCurrentControl(const ModelDataT& data)
        : inductance_(parameter<RealT>(data, InnerCurrentControlParameters::L)),
          kp_(parameter<RealT>(data, InnerCurrentControlParameters::Kp)),
          ki_(parameter<RealT>(data, InnerCurrentControlParameters::Ki)),
          kaw_(parameter<RealT>(data, InnerCurrentControlParameters::Kaw)),
          current_limit_(parameter<RealT>(data, InnerCurrentControlParameters::Imax)),
          modulation_limit_(parameter<RealT>(data, InnerCurrentControlParameters::Mmax)),
          current_coefficient_(RealT{1} / (current_limit_ * current_limit_)),
          voltage_coefficient_(RealT{8} / (RealT{3} * modulation_limit_ * modulation_limit_)),
          monitor_(std::make_unique<MonitorT>(data))
      {
        if (inductance_ <= 0 || kp_ <= 0 || ki_ <= 0 || kaw_ <= 0 || current_limit_ <= 0
            || modulation_limit_ <= 0 || modulation_limit_ > 1
            || current_coefficient_ <= 0 || voltage_coefficient_ <= 0
            || !std::isfinite(current_coefficient_) || !std::isfinite(voltage_coefficient_))
          throw std::invalid_argument("InnerCurrentControl: invalid control parameters");
        this->equation_size_ = this->size_ = 2;
        using Mon                          = typename ModelDataT::MonitorableVariables;
        monitor_->set(Mon::xid, [this]
                      { return this->y_.getData()[0]; });
        monitor_->set(Mon::xiq, [this]
                      { return this->y_.getData()[1]; });
        for (size_t n = 0; n < output_.size(); ++n)
        {
          const auto key = static_cast<Outputs>(n);
          output_[n].setComputed(
              [this, key]
              { return output(key); },
              [this, key](typename SignalT::GradientT& gradient, RealT scale)
              { appendOutputGradient(key, gradient, scale); });
          monitor_->set(static_cast<Mon>(n + 2), [this, key]
                        { return output(key); });
        }
      }

      template <typename scalar_type, typename index_type>
      InnerCurrentControl<scalar_type, index_type>::~InnerCurrentControl() = default;

      template <typename scalar_type, typename index_type>
      typename InnerCurrentControl<scalar_type, index_type>::SignalT&
      InnerCurrentControl<scalar_type, index_type>::outputSignal(Outputs output)
      {
        return output_.at(static_cast<size_t>(output));
      }

      template <typename scalar_type, typename index_type>
      void InnerCurrentControl<scalar_type, index_type>::attachInput(const std::array<SignalT*, 8>& inputs)
      {
        if (this->allocated_)
          throw std::logic_error("InnerCurrentControl: attach inputs before allocation");
        for (auto* signal : inputs)
          if (!signal)
            throw std::invalid_argument("InnerCurrentControl: all inputs are required");
        input_ = inputs;
      }

      template <typename scalar_type, typename index_type>
      void InnerCurrentControl<scalar_type, index_type>::assignOutput(Outputs output, SignalT* signal)
      {
        if (this->allocated_)
          throw std::logic_error("InnerCurrentControl: assign outputs before allocation");
        const auto index    = static_cast<size_t>(output);
        auto&      assigned = alias_.at(index);
        if (!signal || (assigned && assigned != signal))
          throw std::invalid_argument("InnerCurrentControl: invalid output assignment");
        if (assigned == signal)
          return;
        signal->claimProducer();
        assigned = signal;
        signal->setComputed(
            [this, index]
            { return output_[index].read(); },
            [this, index](typename SignalT::GradientT& gradient, RealT scale)
            { output_[index].appendGradient(gradient, scale); });
      }

      template <typename scalar_type, typename index_type>
      int InnerCurrentControl<scalar_type, index_type>::setGridKitComponentID(IdxT id)
      {
        this->gridkit_component_id_ = id;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int InnerCurrentControl<scalar_type, index_type>::allocate()
      {
        for (auto* signal : input_)
          if (!signal)
            throw std::invalid_argument("InnerCurrentControl: all inputs are required");
        if (!this->allocated_)
          this->allocateVectors(2);
        this->tag_.resize(2);
        this->variable_indices_.resize(2);
        this->residual_indices_.resize(2);
        this->assignGlobalIndices(0);
        this->allocateExternalVectors(static_cast<IdxT>(input_.size()), 0);
        for (size_t n = 0; n < input_.size(); ++n)
          this->setExternalVariableSignal(static_cast<IdxT>(n), input_[n]);
        this->allocated_ = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int InnerCurrentControl<scalar_type, index_type>::verify() const
      {
        for (const auto* signal : input_)
          if (!signal || !signal->linked())
            return 1;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int InnerCurrentControl<scalar_type, index_type>::initialize(const std::array<RealT, 2>& integral)
      {
        for (size_t n = 0; n < 2; ++n)
        {
          if (!std::isfinite(integral[n]))
            throw std::invalid_argument("InnerCurrentControl: nonfinite integral state");
          this->y_.getData()[n]  = integral[n];
          this->yp_.getData()[n] = ScalarT{0};
        }
        this->y_.setDataUpdated();
        this->yp_.setDataUpdated();
        return 0;
      }

      template <typename scalar_type, typename index_type>
      void InnerCurrentControl<scalar_type, index_type>::validateInitialState(const std::map<std::string, RealT>& values) const
      {
        for (const auto& [key, value] : values)
          if ((key != "xid" && key != "xiq") || !std::isfinite(value))
            throw std::invalid_argument("InnerCurrentControl: invalid initial state " + key);
      }

      template <typename scalar_type, typename index_type>
      int InnerCurrentControl<scalar_type, index_type>::initializeState(const std::map<std::string, RealT>& values)
      {
        validateInitialState(values);
        std::array<RealT, 2> integral{};
        if (values.contains("xid"))
          integral[0] = values.at("xid");
        if (values.contains("xiq"))
          integral[1] = values.at("xiq");
        return initialize(integral);
      }

      template <typename scalar_type, typename index_type>
      typename InnerCurrentControl<scalar_type, index_type>::Base::InitializationPortsT
      InnerCurrentControl<scalar_type, index_type>::initializationPorts()
      {
        return {};
      }

      template <typename scalar_type, typename index_type>
      int InnerCurrentControl<scalar_type, index_type>::setAbsoluteTolerance(RealT tolerance)
      {
        this->abs_tol_.setToConst(static_cast<ScalarT>(tolerance));
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int InnerCurrentControl<scalar_type, index_type>::evaluateResidual()
      {
        return evaluateInternalResidual();
      }

      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* InnerCurrentControl<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      template <typename scalar_type, typename index_type>
      scalar_type InnerCurrentControl<scalar_type, index_type>::dcVoltage() const
      {
        const auto value = input_[7]->read();
        if (!std::isfinite(static_cast<RealT>(value)) || static_cast<RealT>(value) < RealT{0})
          throw std::invalid_argument("InnerCurrentControl: vdc must be finite and nonnegative");
        return value;
      }

      template <typename scalar_type, typename index_type>
      std::array<scalar_type, 2> InnerCurrentControl<scalar_type, index_type>::limitedReference() const
      {
        const auto d      = input_[4]->read();
        const auto q      = input_[5]->read();
        const auto factor = std::sqrt(Math::max(RealT{1}, current_coefficient_ * (d * d + q * q)));
        return {d / factor, q / factor};
      }

      template <typename scalar_type, typename index_type>
      std::array<scalar_type, 2> InnerCurrentControl<scalar_type, index_type>::unlimitedVoltage() const
      {
        const auto  limited = limitedReference();
        const auto  id      = input_[2]->read();
        const auto  iq      = input_[3]->read();
        const auto  omega_l = input_[6]->read() * inductance_;
        const auto* xi      = this->y_.getData();
        return {input_[0]->read() - omega_l * iq + kp_ * (limited[0] - id) + xi[0],
                input_[1]->read() + omega_l * id + kp_ * (limited[1] - iq) + xi[1]};
      }

      template <typename scalar_type, typename index_type>
      scalar_type InnerCurrentControl<scalar_type, index_type>::output(Outputs output) const
      {
        const auto index = static_cast<size_t>(output);
        if (index >= output_.size() || verify() != 0)
          throw std::logic_error("InnerCurrentControl: invalid output or unconnected input");
        if (index < 2)
          return limitedReference()[index];
        const auto z      = unlimitedVoltage();
        const auto vdc    = dcVoltage();
        const auto factor = std::sqrt(Math::max(vdc * vdc, voltage_coefficient_ * (z[0] * z[0] + z[1] * z[1])));
        return vdc * z.at(index - 2) / factor;
      }

      template <typename scalar_type, typename index_type>
      int InnerCurrentControl<scalar_type, index_type>::evaluateInternalResidual()
      {
        const auto limited = limitedReference();
        const auto z       = unlimitedVoltage();
        for (size_t n = 0; n < 2; ++n)
          this->f_.getData()[n] = ki_ * (limited[n] - input_[n + 2]->read())
                                  + kaw_ * (output(static_cast<Outputs>(n + 2)) - z[n])
                                  - this->yp_.getData()[n];
        this->f_.setDataUpdated();
        return 0;
      }

      template <typename scalar_type, typename index_type>
      void InnerCurrentControl<scalar_type, index_type>::appendLimitedGradient(
          size_t axis, typename SignalT::GradientT& gradient, RealT scale) const
      {
        const std::array<RealT, 2> reference{static_cast<RealT>(input_[4]->read()),
                                             static_cast<RealT>(input_[5]->read())};
        const auto                 norm   = current_coefficient_ * (reference[0] * reference[0] + reference[1] * reference[1]);
        const auto                 factor = std::sqrt(Math::max(RealT{1}, norm));
        const auto                 slope  = Math::sigmoid(norm - RealT{1});
        for (size_t n = 0; n < 2; ++n)
        {
          RealT diagonal = ZERO<RealT>;
          if (n == axis)
            diagonal = ONE<RealT>;
          const auto derivative = diagonal / factor
                                  - reference[axis] * current_coefficient_ * reference[n] * slope
                                        / (factor * factor * factor);
          input_[n + 4]->appendGradient(gradient, scale * derivative);
        }
      }

      template <typename scalar_type, typename index_type>
      void InnerCurrentControl<scalar_type, index_type>::appendUnlimitedVoltageGradient(
          size_t axis, typename SignalT::GradientT& gradient, RealT scale) const
      {
        const auto other = 1 - axis;
        RealT      sign  = ONE<RealT>;
        if (axis == 0)
          sign = -ONE<RealT>;
        input_[axis]->appendGradient(gradient, scale);
        input_[other + 2]->appendGradient(gradient, scale * sign * inductance_ * static_cast<RealT>(input_[6]->read()));
        input_[6]->appendGradient(gradient, scale * sign * inductance_ * static_cast<RealT>(input_[other + 2]->read()));
        appendLimitedGradient(axis, gradient, scale * kp_);
        input_[axis + 2]->appendGradient(gradient, -scale * kp_);
        gradient.emplace_back(this->getVariableIndex(static_cast<IdxT>(axis)), scale);
      }

      template <typename scalar_type, typename index_type>
      void InnerCurrentControl<scalar_type, index_type>::appendOutputGradient(
          Outputs output, typename SignalT::GradientT& gradient, RealT scale) const
      {
        const auto index = static_cast<size_t>(output);
        if (index >= output_.size() || verify() != 0)
          throw std::logic_error("InnerCurrentControl: invalid output or unconnected input");
        if (index < 2)
        {
          appendLimitedGradient(index, gradient, scale);
          return;
        }
        const auto                 axis    = index - 2;
        const auto                 command = unlimitedVoltage();
        const std::array<RealT, 2> z{static_cast<RealT>(command[0]), static_cast<RealT>(command[1])};
        const auto                 vdc          = static_cast<RealT>(dcVoltage());
        const auto                 voltage_norm = voltage_coefficient_ * (z[0] * z[0] + z[1] * z[1]);
        const auto                 factor       = std::sqrt(Math::max(vdc * vdc, voltage_norm));
        const auto                 slope        = Math::sigmoid(voltage_norm - vdc * vdc);
        for (size_t n = 0; n < 2; ++n)
        {
          RealT diagonal = ZERO<RealT>;
          if (n == axis)
            diagonal = ONE<RealT>;
          const auto derivative = vdc * (diagonal / factor - z[axis] * voltage_coefficient_ * z[n] * slope / (factor * factor * factor));
          appendUnlimitedVoltageGradient(n, gradient, scale * derivative);
        }
        const auto derivative = z[axis] / factor
                                - vdc * vdc * z[axis] * (RealT{1} - slope) / (factor * factor * factor);
        input_[7]->appendGradient(gradient, scale * derivative);
      }

      template <typename scalar_type, typename index_type>
      int InnerCurrentControl<scalar_type, index_type>::assembleJacobian(RealT y_scale, RealT yp_scale)
      {
        std::array<typename SignalT::GradientT, 2> gradient;
        size_t                                     entries = 0;
        for (size_t n = 0; n < 2; ++n)
        {
          if (y_scale != RealT{0})
          {
            appendLimitedGradient(n, gradient[n], y_scale * ki_);
            input_[n + 2]->appendGradient(gradient[n], -y_scale * ki_);
            appendOutputGradient(static_cast<Outputs>(n + 2), gradient[n], y_scale * kaw_);
            appendUnlimitedVoltageGradient(n, gradient[n], -y_scale * kaw_);
          }
          if (yp_scale != RealT{0})
            gradient[n].emplace_back(this->getVariableIndex(static_cast<IdxT>(n)), -yp_scale);
          entries += gradient[n].size();
        }
        if (entries != capacity_)
        {
          this->resetJacobianStructure();
          delete[] this->J_rows_buffer_;
          delete[] this->J_cols_buffer_;
          delete[] this->J_vals_buffer_;
          this->J_rows_buffer_ = new IdxT[entries];
          this->J_cols_buffer_ = new IdxT[entries];
          this->J_vals_buffer_ = new RealT[entries];
          capacity_            = entries;
        }
        this->nnz_ = 0;
        for (size_t n = 0; n < 2; ++n)
          for (const auto& [column, value] : gradient[n])
          {
            const auto j            = this->nnz_++;
            this->J_rows_buffer_[j] = this->getResidualIndex(static_cast<IdxT>(n));
            this->J_cols_buffer_[j] = column;
            this->J_vals_buffer_[j] = value;
          }
        if (entries == 0)
          return 0;
        return this->constructCoo();
      }
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
