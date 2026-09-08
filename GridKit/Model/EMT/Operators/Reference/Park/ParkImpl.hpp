#pragma once

#include <cmath>
#include <numbers>
#include <stdexcept>

#include <GridKit/Model/EMT/ComponentInitialization.hpp>
#include <GridKit/Model/EMT/Operators/Reference/Park/Park.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    Park<scalar_type, index_type>::Park()
      : Park(ModelDataT{})
    {
    }

    template <typename scalar_type, typename index_type>
    Park<scalar_type, index_type>::Park(const ModelDataT& data)
      : inverse_(parameter<bool>(data, ModelDataT::Parameters::inverse, false)),
        monitor_(std::make_unique<MonitorT>(data))
    {
      for (size_t n = 0; n < output_port_.size(); ++n)
      {
        const auto key = static_cast<Outputs>(n);
        monitor_->set(static_cast<typename ModelDataT::MonitorableVariables>(n), [this, key]
                      { return output(key); });
        output_port_[n].setComputed(
            [this, key]
            { return output(key); },
            [this, key](typename SignalT::GradientT& gradient, RealT scale)
            { appendOutputGradient(key, gradient, scale); });
      }
    }

    template <typename scalar_type, typename index_type>
    Park<scalar_type, index_type>::~Park() = default;

    template <typename scalar_type, typename index_type>
    int Park<scalar_type, index_type>::setGridKitComponentID(IdxT id)
    {
      this->gridkit_component_id_ = id;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Park<scalar_type, index_type>::allocate()
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
    int Park<scalar_type, index_type>::initialize(const std::map<Outputs, RealT>& outputs)
    {
      if (verify() != 0)
        return 1;
      for (const auto* signal : input_)
      {
        if (!std::isfinite(static_cast<RealT>(signal->read())))
          throw std::invalid_argument("Park inputs must be finite");
      }
      this->validateOutputValues(outputs);
      for (const auto& [key, value] : outputs)
      {
        this->checkOutputValue(outputs, key, static_cast<RealT>(output(key)));
      }
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Park<scalar_type, index_type>::setAbsoluteTolerance(RealT)
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Park<scalar_type, index_type>::evaluateInternalResidual()
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Park<scalar_type, index_type>::evaluateExternalResidual()
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Park<scalar_type, index_type>::evaluateResidual()
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Park<scalar_type, index_type>::assembleJacobian(RealT, RealT)
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* Park<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    template <typename scalar_type, typename index_type>
    void Park<scalar_type, index_type>::attachInput(const std::array<SignalT*, 3>& input, SignalT* theta)
    {
      if (this->allocated_)
        throw std::logic_error("Attach Park inputs before allocation");
      input_ = {input[0], input[1], input[2], theta};
    }

    template <typename scalar_type, typename index_type>
    void Park<scalar_type, index_type>::assignOutput(Outputs output, SignalT* signal)
    {
      const auto index    = static_cast<size_t>(output);
      auto&      assigned = assigned_output_.at(index);
      if (signal == nullptr || (assigned != nullptr && assigned != signal))
        throw std::invalid_argument("Invalid Park output assignment");
      if (assigned == signal)
        return;
      signal->claimProducer();
      assigned = signal;
      signal->setComputed(
          [this, index]
          { return output_port_[index].read(); },
          [this, index](typename SignalT::GradientT& gradient, RealT scale)
          { output_port_[index].appendGradient(gradient, scale); });
    }

    template <typename scalar_type, typename index_type>
    int Park<scalar_type, index_type>::verify() const
    {
      for (const auto* signal : input_)
      {
        if (signal == nullptr || !signal->linked())
        {
          Log::error() << "Park: all vector and angle inputs must have linked sources\n";
          return 1;
        }
      }
      return 0;
    }

    template <typename scalar_type, typename index_type>
    template <typename T>
    ABCMatrix<T> Park<scalar_type, index_type>::transformation(T theta)
    {
      const RealT normalization = std::sqrt(RealT{2} / 3);
      const RealT offset        = 2 * std::numbers::pi_v<RealT> / 3;
      const RealT zero          = 1 / std::sqrt(RealT{3});
      return {{{normalization * std::cos(theta),
                normalization * std::cos(theta - offset),
                normalization * std::cos(theta + offset)},
               {-normalization * std::sin(theta),
                -normalization * std::sin(theta - offset),
                -normalization * std::sin(theta + offset)},
               {T{zero}, T{zero}, T{zero}}}};
    }

    template <typename scalar_type, typename index_type>
    auto Park<scalar_type, index_type>::output(Outputs output) const -> ScalarT
    {
      const auto row = static_cast<size_t>(output);
      if (row >= output_port_.size() || verify() != 0)
        throw std::logic_error("Cannot evaluate an unconnected Park or invalid output");
      const auto matrix = transformation(input_[3]->read());
      ScalarT    value{0};
      for (size_t n = 0; n < 3; ++n)
      {
        auto coefficient = matrix[row][n];
        if (inverse_)
          coefficient = matrix[n][row];
        value += coefficient * input_[n]->read();
      }
      return value;
    }

    template <typename scalar_type, typename index_type>
    void Park<scalar_type, index_type>::appendOutputGradient(
        Outputs output, typename SignalT::GradientT& gradient, RealT scale) const
    {
      const auto row = static_cast<size_t>(output);
      if (row >= output_port_.size() || verify() != 0)
        throw std::logic_error("Cannot differentiate an unconnected Park or invalid output");
      const auto matrix = transformation(static_cast<RealT>(input_[3]->read()));
      RealT      dtheta{0};
      for (size_t n = 0; n < 3; ++n)
      {
        size_t r = row;
        size_t c = n;
        if (inverse_)
        {
          r = n;
          c = row;
        }
        input_[n]->appendGradient(gradient, scale * matrix[r][c]);
        // The angle derivative of the cosine row is the sine row and vice versa.
        RealT derivative = ZERO<RealT>;
        if (r == 0)
          derivative = matrix[1][c];
        if (r == 1)
          derivative = -matrix[0][c];
        dtheta += derivative * static_cast<RealT>(input_[n]->read());
      }
      input_[3]->appendGradient(gradient, scale * dtheta);
    }
  } // namespace EMT
} // namespace GridKit
