#pragma once

#include <limits>
#include <type_traits>
#include <variant>

#include <GridKit/Model/PhasorDynamics/SignalBlock/Delay/Delay.hpp>
#include <GridKit/Model/PhasorDynamics/SignalBlock/Delay/DelayData.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace SignalBlock
    {
      using Log = ::GridKit::Utilities::Logger;

      template <typename scalar_type, typename index_type>
      Delay<scalar_type, index_type>::Delay()
        : monitor_(std::make_unique<MonitorT>())
      {
      }

      template <typename scalar_type, typename index_type>
      Delay<scalar_type, index_type>::Delay(const ModelDataT& data)
        : monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        history_.requireWindow(delay_);
        initializeMonitor();
      }

      template <typename scalar_type, typename index_type>
      Delay<scalar_type, index_type>::~Delay()
      {
      }

      template <typename scalar_type, typename index_type>
      int Delay<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Delay<scalar_type, index_type>::allocate()
      {
        size_ = 0;
        nnz_  = 0;
        y_.clear();
        yp_.clear();
        f_.clear();
        tag_.clear();
        variable_indices_.clear();
        residual_indices_.clear();

        if (signals_.template isAssigned<DelayInternalVariables::OUT>())
        {
          signals_.template getSignalNode<DelayInternalVariables::OUT>()->set(&value_, &invalid_index_);
        }

        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Delay<scalar_type, index_type>::verify() const
      {
        int ret = 0;

        if (!(signals_.template isAttached<DelayExternalVariables::IN>()
              && signals_.template isLinked<DelayExternalVariables::IN>()))
        {
          Log::error() << "Delay: required input signal is not attached/linked\n";
          ret += 1;
        }

        if (!signals_.template isAssigned<DelayInternalVariables::OUT>())
        {
          Log::error() << "Delay: required output signal is not assigned\n";
          ret += 1;
        }
        else if (!signals_.template getSignalNode<DelayInternalVariables::OUT>()->linked())
        {
          Log::error() << "Delay: output signal is assigned but not linked\n";
          ret += 1;
        }

        if (delay_ <= 0.0)
        {
          Log::error() << "Delay: \"delay\" must be positive\n";
          ret += 1;
        }

        return ret;
      }

      template <typename scalar_type, typename index_type>
      int Delay<scalar_type, index_type>::initialize()
      {
        updateTime(time_, 0.0);
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Delay<scalar_type, index_type>::tagDifferentiable()
      {
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Delay<scalar_type, index_type>::evaluateResidual()
      {
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Delay<scalar_type, index_type>::evaluateJacobian()
      {
        return 0;
      }

      template <typename scalar_type, typename index_type>
      void Delay<scalar_type, index_type>::updateTime(RealT t, RealT a)
      {
        Component<ScalarT, IdxT>::updateTime(t, a);
        // Output is a pure function of the past input. The step cap guarantees
        // t - delay_ lies in committed history, so this is explicit forcing.
        value_ = ScalarT{history_.read(t - delay_, prehistory_value_)};
      }

      template <typename scalar_type, typename index_type>
      int Delay<scalar_type, index_type>::resetHistory(RealT t0)
      {
        history_.reset();
        if (signals_.template isAttached<DelayExternalVariables::IN>()
            && signals_.template isLinked<DelayExternalVariables::IN>())
        {
          const RealT u0 = inputValue();
          // Default prehistory holds the input value at the initial condition.
          if (!has_prehistory_)
          {
            prehistory_value_ = u0;
          }
          history_.record(t0, u0);
          value_ = ScalarT{prehistory_value_};
        }
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Delay<scalar_type, index_type>::stepAccepted(RealT t)
      {
        if (signals_.template isAttached<DelayExternalVariables::IN>()
            && signals_.template isLinked<DelayExternalVariables::IN>())
        {
          history_.record(t, inputValue());
        }
        return 0;
      }

      template <typename scalar_type, typename index_type>
      auto Delay<scalar_type, index_type>::maxStepSize() const -> RealT
      {
        return (delay_ > 0.0) ? delay_ : std::numeric_limits<RealT>::infinity();
      }

      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* Delay<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      template <typename scalar_type, typename index_type>
      void Delay<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        auto read_real = [](const auto& value)
        {
          return std::visit([](const auto& stored) -> RealT
                            { return static_cast<RealT>(stored); },
                            value);
        };

        if (data.parameters.contains(ModelDataT::Parameters::delay))
        {
          delay_ = read_real(data.parameters.at(ModelDataT::Parameters::delay));
        }
        if (data.parameters.contains(ModelDataT::Parameters::prehistory))
        {
          prehistory_value_ = read_real(data.parameters.at(ModelDataT::Parameters::prehistory));
          has_prehistory_   = true;
        }
      }

      template <typename scalar_type, typename index_type>
      void Delay<scalar_type, index_type>::initializeMonitor()
      {
        using Variable = typename ModelDataT::MonitorableVariables;
        monitor_->set(Variable::out, [this]
                      { return value_; });
      }

      template <typename scalar_type, typename index_type>
      auto Delay<scalar_type, index_type>::inputValue() const -> RealT
      {
        const ScalarT u = signals_.template readExternalVariable<DelayExternalVariables::IN>();
        if constexpr (std::is_arithmetic_v<ScalarT>)
        {
          return static_cast<RealT>(u);
        }
        else
        {
          return static_cast<RealT>(u.getValue());
        }
      }
    } // namespace SignalBlock
  } // namespace PhasorDynamics
} // namespace GridKit
