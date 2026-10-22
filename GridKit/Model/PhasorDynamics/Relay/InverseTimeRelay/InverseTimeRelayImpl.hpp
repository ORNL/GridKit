/**
 * @file InverseTimeRelayImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the InverseTimeRelay phasor-dynamics relay model.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <format>

#include <GridKit/Model/ConfigurationChecks.hpp>
#include <GridKit/Model/ParameterReader.hpp>
#include <GridKit/Model/PhasorDynamics/Relay/InverseTimeRelay/InverseTimeRelay.hpp>
#include <GridKit/Model/PhasorDynamics/Relay/InverseTimeRelay/InverseTimeRelayData.hpp>
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
       * @brief Construct an InverseTimeRelay from model data
       *
       * @param[in] data Parameters, signal assignments, and monitored-variable
       *                 selections.
       */
      template <typename scalar_type, typename index_type>
      InverseTimeRelay<scalar_type, index_type>::InverseTimeRelay(const ModelDataT& data)
        : monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        initializeMonitor();
        size_ = static_cast<IdxT>(Utilities::enum_size<InverseTimeRelayInternalVariables>());
      }

      /**
       * @brief Destroy the InverseTimeRelay.
       */
      template <typename scalar_type, typename index_type>
      InverseTimeRelay<scalar_type, index_type>::~InverseTimeRelay()
      {
      }

      /**
       * @brief Set the component ID
       *
       * @param[in] component_id Identifier assigned by the system model.
       * @return Zero on success.
       */
      template <typename scalar_type, typename index_type>
      int InverseTimeRelay<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
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
      int InverseTimeRelay<scalar_type, index_type>::allocate()
      {
        const auto S = static_cast<size_t>(InverseTimeRelayInternalVariables::S);

        if (!allocated_)
        {
          this->allocateVectors(size_);
        }
        auto size = static_cast<size_t>(size_);

        tag_.assign(size, false);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        const auto signal_size = Utilities::enum_size<InverseTimeRelayExternalVariables>();
        ws_.resize(static_cast<IdxT>(signal_size));
        ws_.setToZero();
        ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        auto* y = y_.getData();

        if (auto port = ports_.out.template port<InverseTimeRelaySignalOutputs::trip>())
        {
          port.link(&y[S], &(this->getVariableIndex(static_cast<IdxT>(S))));
        }

        allocated_ = true;
        return 0;
      }

      /**
       * @brief Validate the InverseTimeRelay configuration
       *
       * Checks the required parameters, the curve parameters, the required
       * ports, and attached measurement signals.
       *
       * @return The configuration checks; passed() when valid.
       */
      template <typename scalar_type, typename index_type>
      Model::ConfigurationChecks InverseTimeRelay<scalar_type, index_type>::verify() const
      {
        using Inputs = InverseTimeRelaySignalInputs;

        Model::ConfigurationChecks checks("InverseTimeRelay");

        checks.check(has_required_parameters_, "Ipickup, TD, A, B, and alpha are required");
        checks.check(Ipickup_ > ZERO<RealT>, "Ipickup must be positive");
        checks.check(TD_ > ZERO<RealT>, "TD must be positive");
        checks.check(A_ > ZERO<RealT>, "A must be positive");
        checks.check(B_ >= ZERO<RealT>, "B must be non-negative");
        checks.check(exponent_ > ZERO<RealT>, "alpha must be positive");

        checks.check(ports_.in.template port<Inputs::ir>().connected(),
                     "required ir input signal is not attached");
        checks.check(ports_.in.template port<Inputs::ii>().connected(),
                     "required ii input signal is not attached");
        checks.check(ports_.out.template port<InverseTimeRelaySignalOutputs::trip>().connected(),
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

        check_attached_signal.template operator()<Inputs::ir>("ir");
        check_attached_signal.template operator()<Inputs::ii>("ii");

        return checks;
      }

      /**
       * @brief Initialize the relay at reset
       *
       * @return Zero on success; nonzero when the configuration is rejected.
       */
      template <typename scalar_type, typename index_type>
      int InverseTimeRelay<scalar_type, index_type>::initialize()
      {
        const auto X = static_cast<size_t>(InverseTimeRelayInternalVariables::X);
        const auto S = static_cast<size_t>(InverseTimeRelayInternalVariables::S);

        if (!this->verifyAndLog())
        {
          Log::error() << "InverseTimeRelay: cannot initialize with invalid configuration\n";
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
      int InverseTimeRelay<scalar_type, index_type>::tagDifferentiable()
      {
        std::fill(tag_.begin(), tag_.end(), false);
        tag_[static_cast<size_t>(InverseTimeRelayInternalVariables::X)] = true;
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
      int InverseTimeRelay<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      /**
       * @brief Residuals of the relay equations
       *
       * Refreshes the measured-current buffer and evaluates the internal
       * residual.
       *
       * @return Zero on success.
       */
      template <typename scalar_type, typename index_type>
      int InverseTimeRelay<scalar_type, index_type>::evaluateResidual()
      {
        const auto IR = static_cast<size_t>(InverseTimeRelayExternalVariables::IR);
        const auto II = static_cast<size_t>(InverseTimeRelayExternalVariables::II);

        auto* ws = ws_.getData();

        const auto& ir_port = ports_.in.template port<InverseTimeRelaySignalInputs::ir>();
        ws[IR]              = ir_port.readSignal();
        ws_indices_[IR]     = ir_port.signalVariableIndex();

        const auto& ii_port = ports_.in.template port<InverseTimeRelaySignalInputs::ii>();
        ws[II]              = ii_port.readSignal();
        ws_indices_[II]     = ii_port.signalVariableIndex();

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
      const Model::VariableMonitorBase* InverseTimeRelay<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      /**
       * @brief Evaluate the InverseTimeRelay internal residual.
       *
       * The latch rate blends the inverse operate time above pickup with the
       * inverse reset time below it.
       *
       * @param[in] y Internal variables in InverseTimeRelayInternalVariables order.
       * @param[in] yp Internal derivatives in the same order.
       * @param[in] wb Unused; the relay reads no bus.
       * @param[in] ws Measured current in InverseTimeRelayExternalVariables order.
       * @param[out] f Residuals in InverseTimeRelayInternalVariables order.
       * @return Zero on success.
       */
      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline int InverseTimeRelay<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT*                  y,
          const ScalarT*                  yp,
          [[maybe_unused]] const ScalarT* wb,
          const ScalarT*                  ws,
          ScalarT*                        f)
      {
        const auto X = static_cast<size_t>(InverseTimeRelayInternalVariables::X);
        const auto S = static_cast<size_t>(InverseTimeRelayInternalVariables::S);

        const auto IR = static_cast<size_t>(InverseTimeRelayExternalVariables::IR);
        const auto II = static_cast<size_t>(InverseTimeRelayExternalVariables::II);

        const ScalarT x     = y[X];
        const ScalarT s     = y[S];
        const ScalarT x_dot = yp[X];

        const ScalarT ir = ws[IR];
        const ScalarT ii = ws[II];

        const ScalarT m_squared = (ir * ir + ii * ii) / (Ipickup_ * Ipickup_);
        const ScalarT p         = Math::sigmoid(m_squared - ONE<RealT>);

        // K = max(M^2, 1)^(alpha/2) - 1, written with exp and log
        const ScalarT k = std::exp(HALF<RealT> * exponent_ * std::log(Math::max(m_squared, ONE<RealT>)))
                          - ONE<RealT>;

        const ScalarT operate_rate = k / (TD_ * (A_ + B_ * k));
        const ScalarT reset_rate   = (ONE<RealT> - Math::min(m_squared, ONE<RealT>)) / Treset_;
        const ScalarT latch_rate   = LN4 * (p * operate_rate + (ONE<RealT> - p) * reset_rate);

        f[X] = -x_dot + latch_rate * Math::latch(x, p, self_reset_ * (ONE<RealT> - p));
        f[S] = -s + Math::above(x, TRIP_LEVEL);

        return 0;
      }

      /**
       * @brief Read the parameters out of the model data
       *
       * Ipickup, TD, A, B, and alpha are required; a missing one is reported
       * by verify(). A non-numeric or non-finite real, a negative reset
       * constant, or a non-Boolean Lockout throws std::invalid_argument.
       * Integer JSON values are accepted for real parameters.
       *
       * @param[in] data Parameters and monitored-variable selections.
       */
      template <typename scalar_type, typename index_type>
      void InverseTimeRelay<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Params = typename ModelDataT::Parameters;

        const Model::ParameterReader reader(data, "InverseTimeRelay");

        const auto Ipickup = reader.readReal(Params::Ipickup);
        const auto TD      = reader.readReal(Params::TD);
        const auto A       = reader.readReal(Params::A);
        const auto B       = reader.readReal(Params::B);
        const auto alpha   = reader.readReal(Params::alpha);

        has_required_parameters_ = Ipickup.has_value() && TD.has_value() && A.has_value()
                                   && B.has_value() && alpha.has_value();

        Ipickup_  = Ipickup.value_or(Ipickup_);
        TD_       = TD.value_or(TD_);
        A_        = A.value_or(A_);
        B_        = B.value_or(B_);
        exponent_ = alpha.value_or(exponent_);
        tr_       = reader.readNonNegative(Params::tr).value_or(tr_);
        lockout_  = reader.readBool(Params::Lockout).value_or(lockout_);
        setDerivedParameters();
      }

      /**
       * @brief Bind the monitorable variables
       *
       * The measured magnitude is read from the inputs of the last residual
       * evaluation.
       */
      template <typename scalar_type, typename index_type>
      void InverseTimeRelay<scalar_type, index_type>::initializeMonitor()
      {
        using Variable = typename ModelDataT::MonitorableVariables;

        const auto X  = static_cast<size_t>(InverseTimeRelayInternalVariables::X);
        const auto S  = static_cast<size_t>(InverseTimeRelayInternalVariables::S);
        const auto IR = static_cast<size_t>(InverseTimeRelayExternalVariables::IR);
        const auto II = static_cast<size_t>(InverseTimeRelayExternalVariables::II);

        monitor_->set(Variable::im, [this]
                      {
                        const auto* ws = ws_.getData();
                        return std::sqrt(ws[IR] * ws[IR] + ws[II] * ws[II]); });
        monitor_->set(Variable::x, [this]
                      { return y_.getData()[X]; });
        monitor_->set(Variable::trip, [this]
                      { return y_.getData()[S]; });
      }

      /**
       * @brief Resolve the reset time and the lockout mask
       *
       * The time-dial-scaled reset constant is raised to the well-posedness
       * floor, so a zero reset constant resets within the floor.
       */
      template <typename scalar_type, typename index_type>
      void InverseTimeRelay<scalar_type, index_type>::setDerivedParameters()
      {
        Treset_ = std::max(TD_ * tr_, TIME_CONSTANT_MINIMUM);

        self_reset_ = ONE<RealT>;
        if (lockout_)
        {
          self_reset_ = ZERO<RealT>;
        }
      }
    } // namespace Relay
  } // namespace PhasorDynamics
} // namespace GridKit
