/**
 * @file DefiniteTimeRelayImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the DefiniteTimeRelay phasor-dynamics relay model.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <format>

#include <GridKit/Model/ConfigurationChecks.hpp>
#include <GridKit/Model/ParameterReader.hpp>
#include <GridKit/Model/PhasorDynamics/Relay/DefiniteTimeRelay/DefiniteTimeRelay.hpp>
#include <GridKit/Model/PhasorDynamics/Relay/DefiniteTimeRelay/DefiniteTimeRelayData.hpp>
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
      /**
       * @brief Construct a DefiniteTimeRelay from model data
       *
       * @param[in] data Parameters, signal assignments, and monitored-variable
       *                 selections.
       */
      template <typename scalar_type, typename index_type>
      DefiniteTimeRelay<scalar_type, index_type>::DefiniteTimeRelay(const ModelDataT& data)
        : monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        initializeMonitor();
        size_ = static_cast<IdxT>(Utilities::enum_size<DefiniteTimeRelayInternalVariables>());
      }

      /**
       * @brief Destroy the DefiniteTimeRelay.
       */
      template <typename scalar_type, typename index_type>
      DefiniteTimeRelay<scalar_type, index_type>::~DefiniteTimeRelay()
      {
      }

      /**
       * @brief Set the component ID
       *
       * @param[in] component_id Identifier assigned by the system model.
       * @return Zero on success.
       */
      template <typename scalar_type, typename index_type>
      int DefiniteTimeRelay<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      /**
       * @brief Allocate the relay vectors and wire the trip output
       *
       * An assigned `trip` node aliases the trip-command variable.
       *
       * @return Zero on success.
       */
      template <typename scalar_type, typename index_type>
      int DefiniteTimeRelay<scalar_type, index_type>::allocate()
      {
        const auto S = static_cast<size_t>(DefiniteTimeRelayInternalVariables::S);

        if (!allocated_)
        {
          this->allocateVectors(size_);
        }
        auto size = static_cast<size_t>(size_);

        tag_.assign(size, false);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        const auto signal_size = Utilities::enum_size<DefiniteTimeRelayExternalVariables>();
        ws_.resize(static_cast<IdxT>(signal_size));
        ws_.setToZero();
        ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        auto* y = y_.getData();

        if (auto port = ports_.out.template port<DefiniteTimeRelaySignalOutputs::trip>())
        {
          port.link(&y[S], &(this->getVariableIndex(static_cast<IdxT>(S))));
        }

        allocated_ = true;
        return 0;
      }

      /**
       * @brief Validate the DefiniteTimeRelay configuration
       *
       * Checks the required parameters, the pickup and reference levels, the
       * required ports, and attached measurement signals.
       *
       * @return The configuration checks; passed() when valid.
       */
      template <typename scalar_type, typename index_type>
      Model::ConfigurationChecks DefiniteTimeRelay<scalar_type, index_type>::verify() const
      {
        using Inputs = DefiniteTimeRelaySignalInputs;

        Model::ConfigurationChecks checks("DefiniteTimeRelay");

        checks.check(has_required_parameters_, "Mpickup, Mref, and Ttrip are required");
        checks.check(Mpickup_ > ZERO<RealT>, "Mpickup must be positive");
        checks.check(Mref_ >= ZERO<RealT>, "Mref must be non-negative");

        // Equal levels leave no operating direction, so this is an exact
        // comparison by intent rather than a tolerance test.
        checks.check(Mpickup_ != Mref_, "Mpickup must differ from Mref");

        checks.check(ports_.in.template port<Inputs::mr>().connected(),
                     "required mr input signal is not attached");
        checks.check(ports_.out.template port<DefiniteTimeRelaySignalOutputs::trip>().connected(),
                     "required trip output signal is not assigned");

        // An attached port must resolve to writable signal storage. The
        // enumerator is a template argument, so each port names itself once.
        auto check_attached_signal =
            [&]<Inputs variable>(const char* name)
        {
          if (ports_.in.template port<variable>().connected()
              && !ports_.in.template port<variable>().linked())
          {
            checks.fail(std::format("{} signal attached with no linked source", name));
          }
        };

        check_attached_signal.template operator()<Inputs::mr>("mr");
        check_attached_signal.template operator()<Inputs::mi>("mi");

        return checks;
      }

      /**
       * @brief Initialize the relay at reset
       *
       * @return Zero on success; nonzero when the configuration is rejected.
       */
      template <typename scalar_type, typename index_type>
      int DefiniteTimeRelay<scalar_type, index_type>::initialize()
      {
        const auto X = static_cast<size_t>(DefiniteTimeRelayInternalVariables::X);
        const auto S = static_cast<size_t>(DefiniteTimeRelayInternalVariables::S);

        if (!this->verifyAndLog())
        {
          Log::error() << "DefiniteTimeRelay: cannot initialize with invalid configuration\n";
          return 1;
        }

        auto* y = y_.getData();
        y[X]    = ZERO<RealT>;
        y[S]    = Math::above(y[X], TRIP_LEVEL);

        y_.setDataUpdated();
        yp_.setToConst(static_cast<ScalarT>(ZERO<RealT>));

        // For DependencyTracking::Variable, set variable numbers
        if constexpr (std::is_same_v<scalar_type, DependencyTracking::Variable>)
        {
          this->initializeDependencyTrackingVariableNumbers();
        }

        return 0;
      }

      /**
       * @brief Identify the differential variables
       *
       * The relay latch carries a derivative; the trip command is algebraic.
       *
       * @return Zero on success.
       */
      template <typename scalar_type, typename index_type>
      int DefiniteTimeRelay<scalar_type, index_type>::tagDifferentiable()
      {
        std::fill(tag_.begin(), tag_.end(), false);
        tag_[static_cast<size_t>(DefiniteTimeRelayInternalVariables::X)] = true;
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
      int DefiniteTimeRelay<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      /**
       * @brief Residuals of the relay equations
       *
       * Refreshes the measured-quantity buffer and evaluates the internal
       * residual. An unattached `mi` is zero.
       *
       * @return Zero on success.
       */
      template <typename scalar_type, typename index_type>
      int DefiniteTimeRelay<scalar_type, index_type>::evaluateResidual()
      {
        const auto MR = static_cast<size_t>(DefiniteTimeRelayExternalVariables::MR);
        const auto MI = static_cast<size_t>(DefiniteTimeRelayExternalVariables::MI);

        auto* ws = ws_.getData();

        const auto& mr_port = ports_.in.template port<DefiniteTimeRelaySignalInputs::mr>();
        ws[MR]              = mr_port.readSignal();
        ws_indices_[MR]     = mr_port.signalVariableIndex();

        ws[MI]          = ZERO<RealT>;
        ws_indices_[MI] = INVALID_INDEX<IdxT>;
        if (auto mi_port = ports_.in.template port<DefiniteTimeRelaySignalInputs::mi>())
        {
          ws[MI]          = mi_port.readSignal();
          ws_indices_[MI] = mi_port.signalVariableIndex();
        }

        evaluateInternalResidual(y_.getData(), yp_.getData(), nullptr, ws, f_.getData());
        f_.setDataUpdated();
        return 0;
      }

      /**
       * @brief Access the monitor
       *
       * @return Monitor for this model.
       */
      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* DefiniteTimeRelay<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      /**
       * @brief Evaluate the DefiniteTimeRelay internal residual.
       *
       * The pickup indicator normalizes the squared measured magnitude so the
       * reference level maps to -1 and the pickup level to 0 for every element.
       *
       * @param[in] y Internal variables in DefiniteTimeRelayInternalVariables order.
       * @param[in] yp Internal derivatives in the same order.
       * @param[in] wb Unused; the relay reads no bus.
       * @param[in] ws Measured quantity in DefiniteTimeRelayExternalVariables order.
       * @param[out] f Residuals in DefiniteTimeRelayInternalVariables order.
       * @return Zero on success.
       */
      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline int DefiniteTimeRelay<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT*                  y,
          const ScalarT*                  yp,
          [[maybe_unused]] const ScalarT* wb,
          const ScalarT*                  ws,
          ScalarT*                        f)
      {
        const auto X = static_cast<size_t>(DefiniteTimeRelayInternalVariables::X);
        const auto S = static_cast<size_t>(DefiniteTimeRelayInternalVariables::S);

        const auto MR = static_cast<size_t>(DefiniteTimeRelayExternalVariables::MR);
        const auto MI = static_cast<size_t>(DefiniteTimeRelayExternalVariables::MI);

        const ScalarT x     = y[X];
        const ScalarT s     = y[S];
        const ScalarT x_dot = yp[X];

        const ScalarT mr = ws[MR];
        const ScalarT mi = ws[MI];

        const RealT   pickup_squared = Mpickup_ * Mpickup_;
        const ScalarT p              = Math::sigmoid((mr * mr + mi * mi - pickup_squared)
                                        / (pickup_squared - Mref_ * Mref_));

        f[X] = -x_dot + Math::latch(x, p, self_reset_ * (ONE<RealT> - p)) / Tlatch_;
        f[S] = -s + Math::above(x, TRIP_LEVEL);

        return 0;
      }

      /**
       * @brief Read the parameters out of the model data
       *
       * Mpickup, Mref, and Ttrip are required; a missing one is reported by
       * verify(). A non-numeric or non-finite real, a negative trip time, or a
       * non-Boolean Lockout throws std::invalid_argument. Integer JSON values
       * are accepted for real parameters.
       *
       * @param[in] data Parameters and monitored-variable selections.
       */
      template <typename scalar_type, typename index_type>
      void DefiniteTimeRelay<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Params = typename ModelDataT::Parameters;

        const Model::ParameterReader reader(data, "DefiniteTimeRelay");

        const auto Mpickup = reader.readReal(Params::Mpickup);
        const auto Mref    = reader.readReal(Params::Mref);
        const auto Ttrip   = reader.readNonNegative(Params::Ttrip);

        has_required_parameters_ = Mpickup.has_value() && Mref.has_value() && Ttrip.has_value();

        Mpickup_ = Mpickup.value_or(Mpickup_);
        Mref_    = Mref.value_or(Mref_);
        Ttrip_   = Ttrip.value_or(Ttrip_);
        lockout_ = reader.readBool(Params::Lockout).value_or(lockout_);
        setDerivedParameters();
      }

      /**
       * @brief Bind the monitorable variables
       *
       * The measured magnitude is read from the inputs of the last residual
       * evaluation.
       */
      template <typename scalar_type, typename index_type>
      void DefiniteTimeRelay<scalar_type, index_type>::initializeMonitor()
      {
        using Variable = typename ModelDataT::MonitorableVariables;

        const auto X  = static_cast<size_t>(DefiniteTimeRelayInternalVariables::X);
        const auto S  = static_cast<size_t>(DefiniteTimeRelayInternalVariables::S);
        const auto MR = static_cast<size_t>(DefiniteTimeRelayExternalVariables::MR);
        const auto MI = static_cast<size_t>(DefiniteTimeRelayExternalVariables::MI);

        monitor_->set(Variable::m, [this]
                      {
                        const auto* ws = ws_.getData();
                        return std::sqrt(ws[MR] * ws[MR] + ws[MI] * ws[MI]); });
        monitor_->set(Variable::x, [this]
                      { return y_.getData()[X]; });
        monitor_->set(Variable::trip, [this]
                      { return y_.getData()[S]; });
      }

      /**
       * @brief Resolve the latch time and the lockout mask
       *
       * The trip time is raised to the well-posedness floor. The latch time
       * makes the latch reach the trip level after the trip time for a full
       * pickup.
       */
      template <typename scalar_type, typename index_type>
      void DefiniteTimeRelay<scalar_type, index_type>::setDerivedParameters()
      {
        Ttrip_  = std::max(Ttrip_, TIME_CONSTANT_MINIMUM);
        Tlatch_ = Ttrip_ / LN4;

        self_reset_ = ONE<RealT>;
        if (lockout_)
        {
          self_reset_ = ZERO<RealT>;
        }
      }
    } // namespace Relay
  } // namespace PhasorDynamics
} // namespace GridKit
