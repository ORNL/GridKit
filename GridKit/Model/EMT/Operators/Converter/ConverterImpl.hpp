#pragma once

#include <stdexcept>

#include <GridKit/Model/EMT/ComponentInitialization.hpp>
#include <GridKit/Model/EMT/Operators/Converter/Converter.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    Converter<scalar_type, index_type>::Converter()
      : Converter(ModelDataT{})
    {
    }

    template <typename scalar_type, typename index_type>
    Converter<scalar_type, index_type>::Converter(const ModelDataT& data)
      : monitor_(std::make_unique<MonitorT>(data))
    {
      monitor_->set(ModelDataT::MonitorableVariables::voa, [this]
                    { return output(Outputs::voa); });
      monitor_->set(ModelDataT::MonitorableVariables::vob, [this]
                    { return output(Outputs::vob); });
      monitor_->set(ModelDataT::MonitorableVariables::voc, [this]
                    { return output(Outputs::voc); });
      monitor_->set(ModelDataT::MonitorableVariables::idc, [this]
                    { return output(Outputs::idc); });
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
    Converter<scalar_type, index_type>::~Converter() = default;

    template <typename scalar_type, typename index_type>
    void Converter<scalar_type, index_type>::assignOutput(Outputs output, SignalT* signal)
    {
      const auto index    = static_cast<size_t>(output);
      auto&      assigned = assigned_output_.at(index);
      if (signal == nullptr || (assigned != nullptr && assigned != signal))
      {
        throw std::invalid_argument("Invalid Converter output assignment");
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
    int Converter<scalar_type, index_type>::setGridKitComponentID(IdxT id)
    {
      this->gridkit_component_id_ = id;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Converter<scalar_type, index_type>::allocate()
    {
      this->allocateExternalVectors(static_cast<IdxT>(input_.size()), 0);
      for (IdxT n = 0; n < static_cast<IdxT>(input_.size()); ++n)
      {
        this->setExternalVariableSignal(n, input_[static_cast<size_t>(n)]);
      }
      this->allocated_ = true;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Converter<scalar_type, index_type>::initialize(const std::map<Outputs, RealT>& outputs)
    {
      this->validateOutputValues(outputs);
      for (const auto& [key, value] : outputs)
      {
        this->checkOutputValue(outputs, key, static_cast<RealT>(output(key)));
      }
      return verify();
    }

    template <typename scalar_type, typename index_type>
    int Converter<scalar_type, index_type>::setAbsoluteTolerance(RealT)
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Converter<scalar_type, index_type>::evaluateInternalResidual()
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Converter<scalar_type, index_type>::evaluateExternalResidual()
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Converter<scalar_type, index_type>::evaluateResidual()
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Converter<scalar_type, index_type>::assembleJacobian(RealT, RealT)
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* Converter<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    template <typename scalar_type, typename index_type>
    void Converter<scalar_type, index_type>::attachInput(
        const std::array<SignalT*, 3>& switching, SignalT* vdc, const std::array<SignalT*, 3>& current)
    {
      if (this->allocated_)
      {
        throw std::logic_error("Attach Converter inputs before allocation");
      }
      input_ = {switching[0], switching[1], switching[2], vdc, current[0], current[1], current[2]};
    }

    template <typename scalar_type, typename index_type>
    int Converter<scalar_type, index_type>::verify() const
    {
      for (const auto* signal : input_)
      {
        if (signal == nullptr || !signal->linked())
        {
          Log::error() << "Converter: all switching, DC voltage, and AC current inputs must have linked sources\n";
          return 1;
        }
      }
      return 0;
    }

    template <typename scalar_type, typename index_type>
    auto Converter<scalar_type, index_type>::output(Outputs output) const -> ScalarT
    {
      const auto index = static_cast<size_t>(output);
      if (index >= output_port_.size() || verify() != 0)
      {
        throw std::logic_error("Cannot evaluate an unconnected Converter or invalid output");
      }
      const ABCVector<ScalarT> switching{input_[0]->read(), input_[1]->read(), input_[2]->read()};
      if (output == Outputs::idc)
        return dcCurrent(switching, {input_[4]->read(), input_[5]->read(), input_[6]->read()});
      return voltage(switching, input_[3]->read())[index];
    }

    template <typename scalar_type, typename index_type>
    void Converter<scalar_type, index_type>::appendOutputGradient(
        Outputs output, typename SignalT::GradientT& gradient, RealT scale) const
    {
      if (verify() != 0)
      {
        throw std::logic_error("Cannot differentiate an unconnected Converter");
      }
      const auto unit_voltage = voltage({input_[0]->read(), input_[1]->read(), input_[2]->read()}, ScalarT{1});
      if (output == Outputs::idc)
      {
        const auto current = voltage({input_[4]->read(), input_[5]->read(), input_[6]->read()}, ScalarT{1});
        for (size_t n = 0; n < 3; ++n)
        {
          input_[n]->appendGradient(gradient, scale * static_cast<RealT>(current[n]));
          input_[n + 4]->appendGradient(gradient, scale * static_cast<RealT>(unit_voltage[n]));
        }
        return;
      }
      const auto  phase = static_cast<size_t>(output);
      const RealT vdc   = static_cast<RealT>(input_[3]->read());
      for (size_t n = 0; n < 3; ++n)
      {
        input_[n]->appendGradient(gradient, scale * vdc * (n == phase ? RealT{2} : RealT{-1}) / 3);
      }
      input_[3]->appendGradient(gradient, scale * static_cast<RealT>(unit_voltage[phase]));
    }
  } // namespace EMT
} // namespace GridKit
