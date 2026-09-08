#pragma once

#include <cmath>
#include <stdexcept>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/CommonMath.hpp>
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
      initializeParameters(data);
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
    void Modulation<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Parameter = typename ModelDataT::Parameters;
      Mmax_           = parameter<RealT>(data, Parameter::Mmax, Mmax_);
      if (Mmax_ > ZERO<RealT>)
      {
        au_ = RealT{8} / (RealT{3} * Mmax_ * Mmax_);
      }
    }

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
        const std::array<SignalT*, 2>& command, SignalT* vdc)
    {
      if (this->allocated_)
      {
        throw std::logic_error("Attach Modulation inputs before allocation");
      }
      input_ = {command[0], command[1], vdc};
    }

    template <typename scalar_type, typename index_type>
    int Modulation<scalar_type, index_type>::verify() const
    {
      int error_count = 0;
      for (const auto* signal : input_)
      {
        if (signal == nullptr || !signal->linked())
        {
          Log::error() << "Modulation: all voltage command and DC voltage inputs must have linked sources\n";
          ++error_count;
          break;
        }
      }
      if (!std::isfinite(Mmax_) || Mmax_ <= ZERO<RealT> || Mmax_ > ONE<RealT>)
      {
        Log::error() << "Modulation: modulation limit must be finite, positive, and at most one\n";
        ++error_count;
      }
      return error_count;
    }

    template <typename scalar_type, typename index_type>
    auto Modulation<scalar_type, index_type>::output(Outputs output) const -> ScalarT
    {
      const auto index = static_cast<size_t>(output);
      if (index >= output_port_.size() || verify() != 0)
      {
        throw std::logic_error("Cannot evaluate an unconnected Modulation or invalid output");
      }
      const ScalarT vdc = input_[2]->read();
      if (!std::isfinite(static_cast<RealT>(vdc)) || static_cast<RealT>(vdc) < ZERO<RealT>)
        throw std::domain_error("Modulation requires a finite nonnegative DC voltage");
      const ScalarT ud    = input_[0]->read();
      const ScalarT uq    = input_[1]->read();
      const ScalarT scale = std::sqrt(Math::max(vdc * vdc, au_ * (ud * ud + uq * uq)));
      if (output == Outputs::md)
        return ScalarT{2} * ud / scale;
      if (output == Outputs::mq)
        return ScalarT{2} * uq / scale;
      if (output == Outputs::ulimd)
        return vdc * ud / scale;
      return vdc * uq / scale;
    }

    template <typename scalar_type, typename index_type>
    void Modulation<scalar_type, index_type>::appendOutputGradient(
        Outputs output, typename SignalT::GradientT& gradient, RealT scale) const
    {
      if (verify() != 0)
      {
        throw std::logic_error("Cannot differentiate an unconnected Modulation");
      }
      const RealT vdc = static_cast<RealT>(input_[2]->read());
      if (!std::isfinite(vdc) || vdc < ZERO<RealT>)
        throw std::domain_error("Modulation requires a finite nonnegative DC voltage");
      const std::array<RealT, 2> u{static_cast<RealT>(input_[0]->read()), static_cast<RealT>(input_[1]->read())};
      const RealT                radius = au_ * (u[0] * u[0] + u[1] * u[1]);
      const RealT                square = Math::max(vdc * vdc, radius);
      const RealT                norm   = std::sqrt(square);
      // The smooth maximum selects the command radius with a logistic gate.
      const RealT                gate   = Math::sigmoid(radius - vdc * vdc);
      const auto                 index  = static_cast<size_t>(output);
      const auto                 axis   = index % 2;
      RealT                      factor = vdc;
      if (index < 2)
      {
        factor = RealT{2};
      }
      const RealT value = factor * u[axis] / norm;
      for (size_t n = 0; n < 2; ++n)
      {
        RealT partial = -value * gate * au_ * u[n] / square;
        if (n == axis)
        {
          partial += factor / norm;
        }
        input_[n]->appendGradient(gradient, scale * partial);
      }
      RealT dc_partial = -value * (ONE<RealT> - gate) * vdc / square;
      if (index >= 2)
      {
        dc_partial += u[axis] / norm;
      }
      input_[2]->appendGradient(gradient, scale * dc_partial);
    }
  } // namespace EMT
} // namespace GridKit
