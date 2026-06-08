#pragma once

#include <cmath>
#include <type_traits>
#include <variant>

#include <GridKit/Model/PhasorDynamics/SignalBlock/DelaySmooth/DelaySmooth.hpp>
#include <GridKit/Model/PhasorDynamics/SignalBlock/DelaySmooth/DelaySmoothData.hpp>
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
      DelaySmooth<scalar_type, index_type>::DelaySmooth()
        : monitor_(std::make_unique<MonitorT>())
      {
      }

      template <typename scalar_type, typename index_type>
      DelaySmooth<scalar_type, index_type>::DelaySmooth(const ModelDataT& data)
        : monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        initializeMonitor();
      }

      template <typename scalar_type, typename index_type>
      DelaySmooth<scalar_type, index_type>::~DelaySmooth()
      {
      }

      template <typename scalar_type, typename index_type>
      int DelaySmooth<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int DelaySmooth<scalar_type, index_type>::allocate()
      {
        size_     = N_;
        auto size = static_cast<size_t>(size_);

        f_.resize(size);
        y_.resize(size);
        yp_.resize(size);
        tag_.resize(size);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        ws_.resize(1);
        ws_[0] = 0.0;

        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        if (size_ > 0 && signals_.template isAssigned<DelaySmoothInternalVariables::OUT>())
        {
          signals_.template getSignalNode<DelaySmoothInternalVariables::OUT>()->set(
              &y_[size - 1],
              &(this->getVariableIndex(size_ - 1)));
        }

        return 0;
      }

      template <typename scalar_type, typename index_type>
      int DelaySmooth<scalar_type, index_type>::verify() const
      {
        int ret = 0;

        if (!(signals_.template isAttached<DelaySmoothExternalVariables::IN>()
              && signals_.template isLinked<DelaySmoothExternalVariables::IN>()))
        {
          Log::error() << "DelaySmooth: required input signal is not attached/linked\n";
          ret += 1;
        }

        if (!signals_.template isAssigned<DelaySmoothInternalVariables::OUT>())
        {
          Log::error() << "DelaySmooth: required output signal is not assigned\n";
          ret += 1;
        }
        else if (!signals_.template getSignalNode<DelaySmoothInternalVariables::OUT>()->linked())
        {
          Log::error() << "DelaySmooth: output signal is assigned but not linked\n";
          ret += 1;
        }

        if (delay_ <= 0.0)
        {
          Log::error() << "DelaySmooth: \"delay\" must be positive\n";
          ret += 1;
        }
        if (dt_min_ <= 0.0)
        {
          Log::error() << "DelaySmooth: \"dt_min\" must be positive\n";
          ret += 1;
        }
        if (N_ < 1)
        {
          Log::error() << "DelaySmooth: derived block count N must be >= 1\n";
          ret += 1;
        }

        return ret;
      }

      template <typename scalar_type, typename index_type>
      int DelaySmooth<scalar_type, index_type>::initialize()
      {
        RealT u0 = 0.0;
        if (signals_.template isAttached<DelaySmoothExternalVariables::IN>()
            && signals_.template isLinked<DelaySmoothExternalVariables::IN>())
        {
          u0 = inputValue();
        }

        for (IdxT k = 0; k < size_; ++k)
        {
          y_[static_cast<size_t>(k)]  = ScalarT{u0};
          yp_[static_cast<size_t>(k)] = ScalarT{0.0};
        }

        return 0;
      }

      template <typename scalar_type, typename index_type>
      int DelaySmooth<scalar_type, index_type>::tagDifferentiable()
      {
        for (IdxT k = 0; k < size_; ++k)
        {
          tag_[static_cast<size_t>(k)] = true;
        }
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int DelaySmooth<scalar_type, index_type>::evaluateInternalResidual(ScalarT* y, ScalarT* yp, ScalarT* ws, ScalarT* f)
      {
        f[0] = -yp[0] + G_ * (ws[0] - y[0]);
        for (IdxT k = 1; k < N_; ++k)
        {
          f[static_cast<size_t>(k)] =
              -yp[static_cast<size_t>(k)]
              + G_ * (y[static_cast<size_t>(k - 1)] - y[static_cast<size_t>(k)]);
        }
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int DelaySmooth<scalar_type, index_type>::evaluateResidual()
      {
        ws_[0] = 0.0;
        if (signals_.template isAttached<DelaySmoothExternalVariables::IN>())
        {
          ws_[0] = signals_.template readExternalVariable<DelaySmoothExternalVariables::IN>();
        }

        evaluateInternalResidual(y_.data(), yp_.data(), ws_.data(), f_.data());
        return 0;
      }

      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* DelaySmooth<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      template <typename scalar_type, typename index_type>
      void DelaySmooth<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
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
        if (data.parameters.contains(ModelDataT::Parameters::dt_min))
        {
          dt_min_ = read_real(data.parameters.at(ModelDataT::Parameters::dt_min));
        }

        if (delay_ > 0.0 && dt_min_ > 0.0)
        {
          N_ = static_cast<IdxT>(std::floor(delay_ / dt_min_));
          if (N_ >= 1)
          {
            T_ = delay_ / static_cast<RealT>(N_);
            G_ = static_cast<RealT>(1.0) / T_;
          }
        }
      }

      template <typename scalar_type, typename index_type>
      void DelaySmooth<scalar_type, index_type>::initializeMonitor()
      {
        using Variable = typename ModelDataT::MonitorableVariables;
        monitor_->set(Variable::out, [this]
                      { return y_[static_cast<size_t>(N_) - 1]; });
      }

      template <typename scalar_type, typename index_type>
      auto DelaySmooth<scalar_type, index_type>::inputValue() const -> RealT
      {
        const ScalarT u = signals_.template readExternalVariable<DelaySmoothExternalVariables::IN>();
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
