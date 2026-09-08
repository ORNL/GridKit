#pragma once

#include <stdexcept>

#include <GridKit/Model/EMT/ComponentInitialization.hpp>
#include <GridKit/Model/EMT/Operators/Modulation/Modulation.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    Modulation<scalar_type, index_type>::Modulation()
      : Modulation(ModelDataT{})
    {
    }

    template <typename scalar_type, typename index_type>
    Modulation<scalar_type, index_type>::Modulation(const ModelDataT& data)
      : monitor_(std::make_unique<MonitorT>(data))
    {
      for (size_t n = 0; n < output_port_.size(); ++n)
        monitor_->set(static_cast<typename ModelDataT::MonitorableVariables>(n), [this, n]
                      { return output_port_[n].read(); });
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
    Modulation<scalar_type, index_type>::~Modulation() = default;

    template <typename scalar_type, typename index_type>
    void Modulation<scalar_type, index_type>::assignOutput(Outputs output, SignalT* signal)
    {
      const auto index    = static_cast<size_t>(output);
      auto&      assigned = assigned_output_.at(index);
      if (signal == nullptr || (assigned != nullptr && assigned != signal))
      {
        throw std::invalid_argument("Invalid Modulation output assignment");
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
    int Modulation<scalar_type, index_type>::setGridKitComponentID(IdxT id)
    {
      this->gridkit_component_id_ = id;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Modulation<scalar_type, index_type>::allocate()
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
    int Modulation<scalar_type, index_type>::initialize(const std::map<Outputs, RealT>& outputs)
    {
      this->validateOutputValues(outputs);
      for (const auto& [key, value] : outputs)
      {
        this->checkOutputValue(outputs, key, static_cast<RealT>(output(key)));
      }
      return verify();
    }

    template <typename scalar_type, typename index_type>
    int Modulation<scalar_type, index_type>::setAbsoluteTolerance(RealT)
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Modulation<scalar_type, index_type>::evaluateInternalResidual()
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Modulation<scalar_type, index_type>::evaluateExternalResidual()
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Modulation<scalar_type, index_type>::evaluateResidual()
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Modulation<scalar_type, index_type>::assembleJacobian(RealT, RealT)
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* Modulation<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    template <typename scalar_type, typename index_type>
    void Modulation<scalar_type, index_type>::attachInput(
        const std::array<SignalT*, 3>& command, SignalT* vdc)
    {
      if (this->allocated_)
      {
        throw std::logic_error("Attach Modulation inputs before allocation");
      }
      input_ = {command[0], command[1], command[2], vdc};
    }

    template <typename scalar_type, typename index_type>
    int Modulation<scalar_type, index_type>::verify() const
    {
      for (const auto* signal : input_)
      {
        if (signal == nullptr || !signal->linked())
        {
          Log::error() << "Modulation: all voltage command and DC voltage inputs must have linked sources\n";
          return 1;
        }
      }
      return 0;
    }

    template <typename scalar_type, typename index_type>
    auto Modulation<scalar_type, index_type>::output(Outputs output) const -> ScalarT
    {
      const auto index = static_cast<size_t>(output);
      if (index >= output_port_.size() || verify() != 0)
      {
        throw std::logic_error("Cannot evaluate an unconnected Modulation or invalid output");
      }
      const auto vdc = input_[3]->read();
      if (!(static_cast<RealT>(vdc) > RealT{0}) || !std::isfinite(static_cast<RealT>(vdc)))
        throw std::domain_error("Modulation requires a finite positive DC voltage");
      return ScalarT{2} * input_[index]->read() / vdc;
    }

    template <typename scalar_type, typename index_type>
    void Modulation<scalar_type, index_type>::appendOutputGradient(
        Outputs output, typename SignalT::GradientT& gradient, RealT scale) const
    {
      if (verify() != 0)
      {
        throw std::logic_error("Cannot differentiate an unconnected Modulation");
      }
      const auto  phase = static_cast<size_t>(output);
      const RealT vdc   = static_cast<RealT>(input_[3]->read());
      if (!(vdc > RealT{0}) || !std::isfinite(vdc))
        throw std::domain_error("Modulation requires a finite positive DC voltage");
      input_[phase]->appendGradient(gradient, scale * 2 / vdc);
      input_[3]->appendGradient(gradient, -scale * 2 * static_cast<RealT>(input_[phase]->read()) / (vdc * vdc));
    }
  } // namespace EMT
} // namespace GridKit
