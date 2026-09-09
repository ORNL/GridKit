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
      monitor_->set(ModelDataT::MonitorableVariables::ea, [this]
                    { return output_port_[0].read(); });
      monitor_->set(ModelDataT::MonitorableVariables::eb, [this]
                    { return output_port_[1].read(); });
      monitor_->set(ModelDataT::MonitorableVariables::ec, [this]
                    { return output_port_[2].read(); });
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
        const std::array<SignalT*, 3>& switching, SignalT* vdc)
    {
      if (this->allocated_)
      {
        throw std::logic_error("Attach Converter inputs before allocation");
      }
      input_ = {switching[0], switching[1], switching[2], vdc};
    }

    template <typename scalar_type, typename index_type>
    int Converter<scalar_type, index_type>::verify() const
    {
      for (const auto* signal : input_)
      {
        if (signal == nullptr || !signal->linked())
        {
          Log::error() << "Converter: all switching and DC voltage inputs must have linked sources\n";
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
      std::array<ScalarT, 4> values{};
      for (size_t n = 0; n < values.size(); ++n)
        values[n] = input_[n]->read();
      return evaluateOutput(output, values.data());
    }

    template <typename scalar_type, typename index_type>
    auto Converter<scalar_type, index_type>::evaluateOutput(Outputs output, const ScalarT* input) const -> ScalarT
    {
      const ABCVector<ScalarT> switching{input[0], input[1], input[2]};
      return voltage(switching, input[3])[static_cast<size_t>(output)];
    }
  } // namespace EMT
} // namespace GridKit
