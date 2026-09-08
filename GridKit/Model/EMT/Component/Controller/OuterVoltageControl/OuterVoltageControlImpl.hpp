#pragma once

#include <GridKit/Model/EMT/Component/Controller/OuterVoltageControl/OuterVoltageControl.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template <typename scalar_type, typename index_type>
      OuterVoltageControl<scalar_type, index_type>::OuterVoltageControl(const ModelDataT& data)
        : capacitance_(parameter<RealT>(data, OuterVoltageControlParameters::C)),
          kp_(parameter<RealT>(data, OuterVoltageControlParameters::Kp)),
          ki_(parameter<RealT>(data, OuterVoltageControlParameters::Ki)),
          kaw_(parameter<RealT>(data, OuterVoltageControlParameters::Kaw)),
          monitor_(std::make_unique<MonitorT>(data))
      {
        if (capacitance_ <= 0 || kp_ <= 0 || ki_ <= 0 || kaw_ <= 0)
          throw std::invalid_argument("OuterVoltageControl: invalid control parameters");
        this->equation_size_ = this->size_ = 2;
        using Mon                          = typename ModelDataT::MonitorableVariables;
        monitor_->set(Mon::etad, [this]
                      { return this->y_.getData()[0]; });
        monitor_->set(Mon::etaq, [this]
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
      OuterVoltageControl<scalar_type, index_type>::~OuterVoltageControl() = default;

      template <typename scalar_type, typename index_type>
      typename OuterVoltageControl<scalar_type, index_type>::SignalT&
      OuterVoltageControl<scalar_type, index_type>::outputSignal(Outputs output)
      {
        return output_.at(static_cast<size_t>(output));
      }

      template <typename scalar_type, typename index_type>
      void OuterVoltageControl<scalar_type, index_type>::attachInput(const std::array<SignalT*, 9>& inputs)
      {
        if (this->allocated_)
          throw std::logic_error("OuterVoltageControl: attach inputs before allocation");
        for (auto* signal : inputs)
          if (!signal)
            throw std::invalid_argument("OuterVoltageControl: all inputs are required");
        input_ = inputs;
      }

      template <typename scalar_type, typename index_type>
      void OuterVoltageControl<scalar_type, index_type>::assignOutput(Outputs output, SignalT* signal)
      {
        if (this->allocated_)
          throw std::logic_error("OuterVoltageControl: assign outputs before allocation");
        const auto index    = static_cast<size_t>(output);
        auto&      assigned = alias_.at(index);
        if (!signal || (assigned && assigned != signal))
          throw std::invalid_argument("OuterVoltageControl: invalid output assignment");
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
      int OuterVoltageControl<scalar_type, index_type>::setGridKitComponentID(IdxT id)
      {
        this->gridkit_component_id_ = id;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int OuterVoltageControl<scalar_type, index_type>::allocate()
      {
        for (auto* signal : input_)
          if (!signal)
            throw std::invalid_argument("OuterVoltageControl: all inputs are required");
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
      int OuterVoltageControl<scalar_type, index_type>::verify() const
      {
        for (const auto* signal : input_)
          if (!signal || !signal->linked())
            return 1;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int OuterVoltageControl<scalar_type, index_type>::initialize(const std::array<RealT, 2>& integral)
      {
        for (size_t n = 0; n < 2; ++n)
        {
          if (!std::isfinite(integral[n]))
            throw std::invalid_argument("OuterVoltageControl: nonfinite integral state");
          this->y_.getData()[n]  = integral[n];
          this->yp_.getData()[n] = ScalarT{0};
        }
        this->y_.setDataUpdated();
        this->yp_.setDataUpdated();
        return 0;
      }

      template <typename scalar_type, typename index_type>
      void OuterVoltageControl<scalar_type, index_type>::validateInitialState(const std::map<std::string, RealT>& values) const
      {
        for (const auto& [key, value] : values)
          if ((key != "etad" && key != "etaq") || !std::isfinite(value))
            throw std::invalid_argument("OuterVoltageControl: invalid initial state " + key);
      }

      template <typename scalar_type, typename index_type>
      int OuterVoltageControl<scalar_type, index_type>::initializeState(const std::map<std::string, RealT>& values)
      {
        validateInitialState(values);
        std::array<RealT, 2> integral{};
        if (values.contains("etad"))
          integral[0] = values.at("etad");
        if (values.contains("etaq"))
          integral[1] = values.at("etaq");
        return initialize(integral);
      }

      template <typename scalar_type, typename index_type>
      typename OuterVoltageControl<scalar_type, index_type>::Base::InitializationPortsT
      OuterVoltageControl<scalar_type, index_type>::initializationPorts()
      {
        return {};
      }

      template <typename scalar_type, typename index_type>
      int OuterVoltageControl<scalar_type, index_type>::setAbsoluteTolerance(RealT tolerance)
      {
        this->abs_tol_.setToConst(static_cast<ScalarT>(tolerance));
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int OuterVoltageControl<scalar_type, index_type>::evaluateResidual()
      {
        return evaluateInternalResidual();
      }

      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* OuterVoltageControl<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      template <typename scalar_type, typename index_type>
      scalar_type OuterVoltageControl<scalar_type, index_type>::output(Outputs output) const
      {
        const auto axis = static_cast<size_t>(output);
        if (axis >= output_.size() || verify() != 0)
          throw std::logic_error("OuterVoltageControl: invalid output or unconnected input");
        const auto other = 1 - axis;
        RealT      sign  = ONE<RealT>;
        if (axis == 0)
          sign = -ONE<RealT>;
        return input_[axis + 4]->read()
               + sign * input_[6]->read() * capacitance_ * input_[other + 2]->read()
               + kp_ * (input_[axis]->read() - input_[axis + 2]->read())
               + this->y_.getData()[axis];
      }

      template <typename scalar_type, typename index_type>
      int OuterVoltageControl<scalar_type, index_type>::evaluateInternalResidual()
      {
        for (size_t n = 0; n < 2; ++n)
          this->f_.getData()[n] = ki_ * (input_[n]->read() - input_[n + 2]->read())
                                  + kaw_ * (input_[n + 7]->read() - output(static_cast<Outputs>(n)))
                                  - this->yp_.getData()[n];
        this->f_.setDataUpdated();
        return 0;
      }

      template <typename scalar_type, typename index_type>
      void OuterVoltageControl<scalar_type, index_type>::appendOutputGradient(
          Outputs output, typename SignalT::GradientT& gradient, RealT scale) const
      {
        const auto axis = static_cast<size_t>(output);
        if (axis >= output_.size() || verify() != 0)
          throw std::logic_error("OuterVoltageControl: invalid output or unconnected input");
        const auto other = 1 - axis;
        RealT      sign  = ONE<RealT>;
        if (axis == 0)
          sign = -ONE<RealT>;
        input_[axis + 4]->appendGradient(gradient, scale);
        input_[other + 2]->appendGradient(gradient, scale * sign * capacitance_ * static_cast<RealT>(input_[6]->read()));
        input_[6]->appendGradient(gradient, scale * sign * capacitance_ * static_cast<RealT>(input_[other + 2]->read()));
        input_[axis]->appendGradient(gradient, scale * kp_);
        input_[axis + 2]->appendGradient(gradient, -scale * kp_);
        gradient.emplace_back(this->getVariableIndex(static_cast<IdxT>(axis)), scale);
      }

      template <typename scalar_type, typename index_type>
      int OuterVoltageControl<scalar_type, index_type>::assembleJacobian(RealT y_scale, RealT yp_scale)
      {
        std::array<typename SignalT::GradientT, 2> gradient;
        size_t                                     entries = 0;
        for (size_t n = 0; n < 2; ++n)
        {
          if (y_scale != RealT{0})
          {
            input_[n]->appendGradient(gradient[n], y_scale * ki_);
            input_[n + 2]->appendGradient(gradient[n], -y_scale * ki_);
            input_[n + 7]->appendGradient(gradient[n], y_scale * kaw_);
            appendOutputGradient(static_cast<Outputs>(n), gradient[n], -y_scale * kaw_);
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
