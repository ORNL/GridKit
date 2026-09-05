#pragma once

#include <stdexcept>

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
                    { return output(0); });
      monitor_->set(ModelDataT::MonitorableVariables::vob, [this]
                    { return output(1); });
      monitor_->set(ModelDataT::MonitorableVariables::voc, [this]
                    { return output(2); });
      for (size_t phase = 0; phase < 3; ++phase)
      {
        output_port_.signals[phase].setComputed(
            [this, phase]
            { return output(phase); },
            [this, phase](typename SignalT::GradientT& gradient, RealT scale)
            { appendOutputGradient(phase, gradient, scale); });
      }
    }

    template <typename scalar_type, typename index_type>
    Converter<scalar_type, index_type>::~Converter() = default;

    template <typename scalar_type, typename index_type>
    void Converter<scalar_type, index_type>::assignOutput(size_t phase, SignalT* signal)
    {
      auto& assigned = assigned_output_.at(phase);
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
          [this, phase]
          { return output_port_.signals[phase].read(); },
          [this, phase](typename SignalT::GradientT& gradient, RealT scale)
          { output_port_.signals[phase].appendGradient(gradient, scale); });
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
      this->allocateExternalVectors(4, 0);
      for (IdxT n = 0; n < 4; ++n)
      {
        this->setExternalVariableSignal(n, input_[static_cast<size_t>(n)]);
      }
      this->allocated_ = true;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Converter<scalar_type, index_type>::initialize()
    {
      return verify();
    }

    template <typename scalar_type, typename index_type>
    int Converter<scalar_type, index_type>::tagDifferentiable()
    {
      return 0;
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
    int Converter<scalar_type, index_type>::evaluateJacobian()
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* Converter<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    template <typename scalar_type, typename index_type>
    void Converter<scalar_type, index_type>::attachInput(SignalT* a, SignalT* b, SignalT* c, SignalT* vdc)
    {
      if (this->allocated_)
      {
        throw std::logic_error("Attach Converter inputs before allocation");
      }
      input_ = {a, b, c, vdc};
    }

    template <typename scalar_type, typename index_type>
    void Converter<scalar_type, index_type>::attachInput(Port3T* s, SignalT* vdc)
    {
      if (s == nullptr)
      {
        throw std::invalid_argument("Converter switching port is null");
      }
      attachInput(s->a(), s->b(), s->c(), vdc);
    }

    template <typename scalar_type, typename index_type>
    int Converter<scalar_type, index_type>::verify() const
    {
      for (const auto* signal : input_)
      {
        if (signal == nullptr || !signal->linked())
        {
          Log::error() << "Converter: all switching and DC-link inputs must have linked sources\n";
          return 1;
        }
      }
      return 0;
    }

    template <typename scalar_type, typename index_type>
    auto Converter<scalar_type, index_type>::output(size_t phase) const -> ScalarT
    {
      if (phase >= 3 || verify() != 0)
      {
        throw std::logic_error("Cannot evaluate an unconnected Converter or invalid phase");
      }
      return voltage({input_[0]->read(), input_[1]->read(), input_[2]->read()}, input_[3]->read())[phase];
    }

    template <typename scalar_type, typename index_type>
    void Converter<scalar_type, index_type>::appendOutputGradient(
        size_t phase, typename SignalT::GradientT& gradient, RealT scale) const
    {
      if (verify() != 0)
      {
        throw std::logic_error("Cannot differentiate an unconnected Converter");
      }
      const RealT vdc = static_cast<RealT>(input_[3]->read());
      for (size_t n = 0; n < 3; ++n)
      {
        input_[n]->appendGradient(gradient, scale * vdc * (n == phase ? RealT{2} : RealT{-1}) / 3);
      }
      const auto unit_voltage = voltage({input_[0]->read(), input_[1]->read(), input_[2]->read()}, ScalarT{1});
      input_[3]->appendGradient(gradient, scale * static_cast<RealT>(unit_voltage[phase]));
    }
  } // namespace EMT
} // namespace GridKit
