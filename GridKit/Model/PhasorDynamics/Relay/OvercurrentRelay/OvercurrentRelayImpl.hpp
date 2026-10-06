/**
 * @file OvercurrentRelayImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the overcurrent relay model.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <variant>

#include <GridKit/Model/PhasorDynamics/Relay/OvercurrentRelay/OvercurrentRelay.hpp>
#include <GridKit/Model/PhasorDynamics/Relay/OvercurrentRelay/OvercurrentRelayData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Enum.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Relay
    {
      /// Logger used for overcurrent relay diagnostics.
      using Log = ::GridKit::Utilities::Logger;

      /**
       * @brief Construct an overcurrent relay without parameters
       *
       * The model is sized but left unconfigured, and no monitor is created.
       */
      template <typename scalar_type, typename index_type>
      OvercurrentRelay<scalar_type, index_type>::OvercurrentRelay()
      {
        size_ = static_cast<IdxT>(Utilities::enum_size<OvercurrentRelayInternalVariables>());
        setDerivedParameters();
      }

      /**
       * @brief Construct an overcurrent relay from model data
       *
       * @param[in] data Parameters, ports, and monitored-variable selections.
       */
      template <typename scalar_type, typename index_type>
      OvercurrentRelay<scalar_type, index_type>::OvercurrentRelay(const ModelDataT& data)
        : monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        initializeMonitor();
        size_ = static_cast<IdxT>(Utilities::enum_size<OvercurrentRelayInternalVariables>());
      }

      /**
       * @brief Destroy the overcurrent relay.
       */
      template <typename scalar_type, typename index_type>
      OvercurrentRelay<scalar_type, index_type>::~OvercurrentRelay()
      {
      }

      /**
       * @brief Set the component ID
       *
       * @param[in] component_id Identifier assigned by the system model.
       * @return Zero on success.
       */
      template <typename scalar_type, typename index_type>
      int OvercurrentRelay<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      /**
       * @brief Allocate the model vectors and wire the trip output
       *
       * @return Zero on success.
       */
      template <typename scalar_type, typename index_type>
      int OvercurrentRelay<scalar_type, index_type>::allocate()
      {
        const auto TRIP = static_cast<size_t>(OvercurrentRelayInternalVariables::TRIP);

        if (!allocated_)
        {
          this->allocateVectors(size_);
        }
        auto size = static_cast<size_t>(size_);

        tag_.assign(size, false);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        const auto signal_size = Utilities::enum_size<OvercurrentRelayExternalVariables>();
        ws_.resize(static_cast<IdxT>(signal_size));
        ws_.setToZero();
        ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        if (auto port = ports_.out.template port<OvercurrentRelaySignalOutputs::trip>())
        {
          port.link(&y_.getData()[TRIP], &(this->getVariableIndex(static_cast<IdxT>(TRIP))));
        }

        allocated_ = true;
        return 0;
      }

      /**
       * @brief Validate the overcurrent relay configuration
       *
       * @return Number of configuration errors; zero when valid.
       */
      template <typename scalar_type, typename index_type>
      int OvercurrentRelay<scalar_type, index_type>::verify() const
      {
        int ret = static_cast<int>(parameter_error_count_);

        auto check = [&](bool condition, const char* message)
        {
          if (!condition)
          {
            Log::error() << "OvercurrentRelay: " << message << '\n';
            ret += 1;
          }
        };

        check(Ipickup_ > ZERO<RealT>, "Ipickup must be positive");
        check(ports_.in.template port<OvercurrentRelaySignalInputs::ir>().linked(),
              "required ir input signal is not linked");
        check(ports_.in.template port<OvercurrentRelaySignalInputs::ii>().linked(),
              "required ii input signal is not linked");
        check(ports_.out.template port<OvercurrentRelaySignalOutputs::trip>().connected(),
              "required trip output signal is not assigned");

        return ret;
      }

      /**
       * @brief Initialize the relay reset
       *
       * @return Zero on success; nonzero when the configuration is rejected.
       */
      template <typename scalar_type, typename index_type>
      int OvercurrentRelay<scalar_type, index_type>::initialize()
      {
        const auto X    = static_cast<size_t>(OvercurrentRelayInternalVariables::X);
        const auto TRIP = static_cast<size_t>(OvercurrentRelayInternalVariables::TRIP);

        if (!allocated_)
        {
          Log::error() << "OvercurrentRelay: allocate must complete before initialize\n";
          return 1;
        }

        if (verify() > 0)
        {
          Log::error() << "OvercurrentRelay: cannot initialize with invalid configuration\n";
          return 1;
        }

        auto* y = y_.getData();
        y[X]    = ZERO<RealT>;
        y[TRIP] = Math::above(y[X], TRIP_LEVEL);

        y_.setDataUpdated();
        yp_.setToConst(static_cast<ScalarT>(ZERO<RealT>));

        if constexpr (std::is_same_v<scalar_type, DependencyTracking::Variable>)
        {
          this->initializeDependencyTrackingVariableNumbers();
        }

        return 0;
      }

      /**
       * @brief Identify the differential variables
       *
       * The latch carries a derivative; the trip command is algebraic.
       *
       * @return Zero on success.
       */
      template <typename scalar_type, typename index_type>
      int OvercurrentRelay<scalar_type, index_type>::tagDifferentiable()
      {
        std::fill(tag_.begin(), tag_.end(), false);
        tag_[static_cast<size_t>(OvercurrentRelayInternalVariables::X)] = true;
        return 0;
      }

      /**
       * @brief Compute the absolute tolerance for each variable in the model
       *
       * Every internal variable receives @p rel_tol as its absolute tolerance.
       *
       * @param[in] rel_tol Solver relative tolerance.
       * @return Zero on success.
       */
      template <typename scalar_type, typename index_type>
      int OvercurrentRelay<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      /**
       * @brief Residuals of system equations
       *
       * Refreshes the signal interface buffer and evaluates the internal residual.
       *
       * @return Zero on success.
       */
      template <typename scalar_type, typename index_type>
      int OvercurrentRelay<scalar_type, index_type>::evaluateResidual()
      {
        const auto IR = static_cast<size_t>(OvercurrentRelayExternalVariables::IR);
        const auto II = static_cast<size_t>(OvercurrentRelayExternalVariables::II);

        auto* ws = ws_.getData();

        ws[IR]          = ports_.in.template port<OvercurrentRelaySignalInputs::ir>().readSignal();
        ws_indices_[IR] = ports_.in.template port<OvercurrentRelaySignalInputs::ir>().signalVariableIndex();
        ws[II]          = ports_.in.template port<OvercurrentRelaySignalInputs::ii>().readSignal();
        ws_indices_[II] = ports_.in.template port<OvercurrentRelaySignalInputs::ii>().signalVariableIndex();

        evaluateInternalResidual(y_.getData(), yp_.getData(), nullptr, ws, f_.getData());
        f_.setDataUpdated();
        return 0;
      }

      /**
       * @brief Access the monitor
       *
       * @return Monitor for this model, or nullptr when the model was
       *         constructed without data.
       */
      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* OvercurrentRelay<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      /**
       * @brief Evaluate the lockout latch and the trip command
       *
       * The body is kept free of branches and loops so sparse automatic
       * differentiation resolves a fixed structure.
       *
       * @param[in] y Internal variables in OvercurrentRelayInternalVariables order.
       * @param[in] yp Internal derivatives in the same enum order.
       * @param[in] wb Unused; the relay has no bus.
       * @param[in] ws Signal values in OvercurrentRelayExternalVariables order.
       * @param[out] f Residuals in OvercurrentRelayInternalVariables order.
       * @return Zero on success.
       */
      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline int
      OvercurrentRelay<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT*                  y,
          const ScalarT*                  yp,
          [[maybe_unused]] const ScalarT* wb,
          const ScalarT*                  ws,
          ScalarT*                        f)
      {
        const auto X    = static_cast<size_t>(OvercurrentRelayInternalVariables::X);
        const auto TRIP = static_cast<size_t>(OvercurrentRelayInternalVariables::TRIP);

        const auto IR = static_cast<size_t>(OvercurrentRelayExternalVariables::IR);
        const auto II = static_cast<size_t>(OvercurrentRelayExternalVariables::II);

        const ScalarT x    = y[X];
        const ScalarT trip = y[TRIP];

        const ScalarT x_dot = yp[X];

        const ScalarT ir = ws[IR];
        const ScalarT ii = ws[II];

        const ScalarT pickup = Math::above((ir * ir + ii * ii) * inv_Ipickup2_, ONE<RealT>);

        f[X]    = -x_dot + Math::latch(x, pickup, ZERO<RealT>) / Tlatch_;
        f[TRIP] = -trip + Math::above(x, TRIP_LEVEL);

        return 0;
      }

      //
      //  Private methods
      //

      /**
       * @brief Read the parameters out of the model data
       *
       * Both settings are required. A missing, non-numeric, or non-finite
       * value is counted and reported by verify() rather than throwing.
       * Integer JSON values are accepted for real parameters.
       *
       * @param[in] data Parameters and monitored-variable selections.
       */
      template <typename scalar_type, typename index_type>
      void OvercurrentRelay<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Params = typename ModelDataT::Parameters;

        parameter_error_count_ = 0;

        auto load_real = [&](auto key, RealT& target, const char* name)
        {
          if (!data.parameters.contains(key))
          {
            Log::error() << "OvercurrentRelay: parameter '" << name << "' is required\n";
            ++parameter_error_count_;
            return;
          }

          const auto& value = data.parameters.at(key);
          RealT       parsed_value{};
          if (const auto* real_value = std::get_if<RealT>(&value))
          {
            parsed_value = *real_value;
          }
          else if (const auto* index_value = std::get_if<IdxT>(&value))
          {
            parsed_value = static_cast<RealT>(*index_value);
          }
          else
          {
            Log::error() << "OvercurrentRelay: parameter '" << name << "' must be numeric\n";
            ++parameter_error_count_;
            return;
          }

          if (!std::isfinite(parsed_value))
          {
            Log::error() << "OvercurrentRelay: parameter '" << name << "' must be finite\n";
            ++parameter_error_count_;
            return;
          }

          target = parsed_value;
        };

        load_real(Params::Ipickup, Ipickup_, "Ipickup");
        load_real(Params::Ttrip, Ttrip_, "Ttrip");
        setDerivedParameters();
      }

      /**
       * @brief Bind the monitorable variables
       */
      template <typename scalar_type, typename index_type>
      void OvercurrentRelay<scalar_type, index_type>::initializeMonitor()
      {
        using Variable = typename ModelDataT::MonitorableVariables;

        monitor_->set(Variable::im, [this]
                      {
                        const ScalarT ir = ports_.in.template port<OvercurrentRelaySignalInputs::ir>().readSignal();
                        const ScalarT ii = ports_.in.template port<OvercurrentRelaySignalInputs::ii>().readSignal();
                        return std::sqrt(ir * ir + ii * ii); });
        monitor_->set(Variable::x, [this]
                      { return y_.getData()[static_cast<size_t>(OvercurrentRelayInternalVariables::X)]; });
        monitor_->set(Variable::trip, [this]
                      { return y_.getData()[static_cast<size_t>(OvercurrentRelayInternalVariables::TRIP)]; });
      }

      /**
       * @brief Static method to log time constant warnings
       */
      template <typename scalar_type, typename index_type>
      void OvercurrentRelay<scalar_type, index_type>::logTimeConstantWarning()
      {
        Log::warning() << "OvercurrentRelay: Ttrip below "
                       << TIME_CONSTANT_MINIMUM
                       << " s is raised to that floor to keep the latch well posed\n";
      }

      /**
       * @brief Resolve the pickup scaling and the latch time constant
       */
      template <typename scalar_type, typename index_type>
      void OvercurrentRelay<scalar_type, index_type>::setDerivedParameters()
      {
        if (Ttrip_ < ZERO<RealT>)
        {
          Log::error() << "OvercurrentRelay: Ttrip must be non-negative\n";
          ++parameter_error_count_;
        }
        if (Ttrip_ < TIME_CONSTANT_MINIMUM)
        {
          logTimeConstantWarning();
          Ttrip_ = TIME_CONSTANT_MINIMUM;
        }

        inv_Ipickup2_ = ONE<RealT> / (Ipickup_ * Ipickup_);

        // Full pickup drives x = 1 - exp(-t / Tlatch) to TRIP_LEVEL in Ttrip.
        Tlatch_ = -Ttrip_ / std::log1p(-TRIP_LEVEL);
      }
    } // namespace Relay
  } // namespace PhasorDynamics
} // namespace GridKit
