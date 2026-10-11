/**
 * @file IeeestImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the IEEEST Power System Stabilizer.
 */

#pragma once

#include <algorithm>
#include <atomic>
#include <format>

#include <GridKit/Model/ConfigurationChecks.hpp>
#include <GridKit/Model/ParameterReader.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/Ieeest.hpp>
#include <GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/IeeestData.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Enum.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Stabilizer
    {
      /// Logger used for IEEEST diagnostics.
      using Log = ::GridKit::Utilities::Logger;

      /**
       * @brief Construct an IEEEST stabilizer without parameters
       *
       * Every parameter keeps its documented default, and no monitor is
       * created. verify() rejects the model until an input is attached.
       */
      template <typename scalar_type, typename index_type, size_t order>
      Ieeest<scalar_type, index_type, order>::Ieeest()
      {
        size_ = static_cast<IdxT>(static_cast<size_t>(IeeestInternalVariables::X1) + order);
        setDerivedParameters();
      }

      /**
       * @brief Construct an IEEEST stabilizer from model data
       *
       * @param[in] data Parameters and monitored-variable selections.
       */
      template <typename scalar_type, typename index_type, size_t order>
      Ieeest<scalar_type, index_type, order>::Ieeest(const ModelDataT& data)
        : monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        initializeMonitor();
        size_ = static_cast<IdxT>(static_cast<size_t>(IeeestInternalVariables::X1) + order);
      }

      /**
       * @brief Destroy the IEEEST stabilizer.
       */
      template <typename scalar_type, typename index_type, size_t order>
      Ieeest<scalar_type, index_type, order>::~Ieeest()
      {
      }

      /**
       * @brief Set the component ID
       *
       * @param[in] component_id Identifier assigned by the system model.
       * @return Zero on success.
       */
      template <typename scalar_type, typename index_type, size_t order>
      int Ieeest<scalar_type, index_type, order>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      /**
       * @brief Allocate the model vectors and wire the stabilizer output
       *
       * An assigned `output` node aliases the limited stabilizer signal.
       * Repeated calls reuse the allocated vectors.
       *
       * @return Zero on success.
       */
      template <typename scalar_type, typename index_type, size_t order>
      int Ieeest<scalar_type, index_type, order>::allocate()
      {
        const auto VSS = static_cast<size_t>(IeeestInternalVariables::VSS);

        if (!allocated_)
        {
          this->allocateVectors(size_);
        }
        const auto size = static_cast<size_t>(size_);

        tag_.assign(size, false);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        const auto signal_size = Utilities::enum_size<IeeestExternalVariables>();
        ws_.resize(static_cast<IdxT>(signal_size));
        ws_.setToZero();
        ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        if (auto port = ports_.out.template port<IeeestSignalOutputs::output>())
        {
          auto* y = y_.getData();
          port.link(&y[VSS], &(this->getVariableIndex(static_cast<IdxT>(VSS))));
        }

        allocated_ = true;
        return 0;
      }

      /**
       * @brief Validate the IEEEST configuration
       *
       * Checks the input signal, the output limits, and that the notch
       * coefficients define this specialization's order with a numerator of
       * no higher degree.
       *
       * @return The configuration checks; passed() when valid.
       */
      template <typename scalar_type, typename index_type, size_t order>
      Model::ConfigurationChecks Ieeest<scalar_type, index_type, order>::verify() const
      {
        Model::ConfigurationChecks checks("Ieeest");

        const auto input = ports_.in.template port<IeeestSignalInputs::input>();
        if (!input.connected())
        {
          checks.fail("required input signal is not attached");
        }
        if (input.connected() && !input.linked())
        {
          checks.fail("input signal attached with no linked source");
        }

        checks.check(Lsmin_ < Lsmax_, "Lsmin must be less than Lsmax");

        // Model data uses an exact zero to mark an absent notch term, so the
        // degrees come from exact comparisons by intent rather than tolerance tests.
        const size_t data_order = notchOrder(A1_, A2_, A3_, A4_);
        if (data_order != order)
        {
          checks.fail(std::format("A1 to A4 define notch order {}, not {}", data_order, order));
        }
        checks.check(factorDegree(A5_, A6_) <= order,
                     "notch numerator degree must not exceed the notch order");

        return checks;
      }

      /**
       * @brief Initialize IEEEST at the steady state of its input
       *
       * A constant input is an equilibrium: each notch factor and lead-lag
       * passes it through, and the washout blocks it, so the output is zero.
       *
       * @return Zero on success; nonzero when the configuration is rejected.
       */
      template <typename scalar_type, typename index_type, size_t order>
      int Ieeest<scalar_type, index_type, order>::initialize()
      {
        const auto X5  = static_cast<size_t>(IeeestInternalVariables::X5);
        const auto X6  = static_cast<size_t>(IeeestInternalVariables::X6);
        const auto X7  = static_cast<size_t>(IeeestInternalVariables::X7);
        const auto VSS = static_cast<size_t>(IeeestInternalVariables::VSS);
        const auto X1  = static_cast<size_t>(IeeestInternalVariables::X1);
        const auto X2  = static_cast<size_t>(IeeestInternalVariables::X2);
        const auto X3  = static_cast<size_t>(IeeestInternalVariables::X3);
        const auto X4  = static_cast<size_t>(IeeestInternalVariables::X4);
        const auto U   = static_cast<size_t>(IeeestExternalVariables::U);

        if (!this->verifyAndLog())
        {
          Log::error() << "Ieeest: cannot initialize with invalid configuration\n";
          return 1;
        }

        auto          input = ports_.in.template port<IeeestSignalInputs::input>();
        const ScalarT u     = input.readSignal();

        ws_.getData()[U] = u;
        ws_indices_[U]   = input.signalVariableIndex();

        auto* y = y_.getData();

        y[X5]  = u;
        y[X6]  = u;
        y[X7]  = u;
        y[VSS] = Math::clamp(ZERO<RealT>, Lsmin_, Lsmax_);
        if constexpr (order >= 1)
        {
          y[X1] = u;
        }
        if constexpr (order >= 2)
        {
          y[X2] = ZERO<RealT>;
        }
        if constexpr (order >= 3)
        {
          y[X3] = u;
        }
        if constexpr (order == 4)
        {
          y[X4] = ZERO<RealT>;
        }

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
       * Every internal variable is differential except the limited output.
       *
       * @return Zero on success.
       */
      template <typename scalar_type, typename index_type, size_t order>
      int Ieeest<scalar_type, index_type, order>::tagDifferentiable()
      {
        const auto X5  = static_cast<size_t>(IeeestInternalVariables::X5);
        const auto X6  = static_cast<size_t>(IeeestInternalVariables::X6);
        const auto X7  = static_cast<size_t>(IeeestInternalVariables::X7);
        const auto VSS = static_cast<size_t>(IeeestInternalVariables::VSS);
        const auto X1  = static_cast<size_t>(IeeestInternalVariables::X1);
        const auto X2  = static_cast<size_t>(IeeestInternalVariables::X2);
        const auto X3  = static_cast<size_t>(IeeestInternalVariables::X3);
        const auto X4  = static_cast<size_t>(IeeestInternalVariables::X4);

        tag_[X5]  = true;
        tag_[X6]  = true;
        tag_[X7]  = true;
        tag_[VSS] = false;
        if constexpr (order >= 1)
        {
          tag_[X1] = true;
        }
        if constexpr (order >= 2)
        {
          tag_[X2] = true;
        }
        if constexpr (order >= 3)
        {
          tag_[X3] = true;
        }
        if constexpr (order == 4)
        {
          tag_[X4] = true;
        }
        return 0;
      }

      /**
       * @brief Compute the absolute tolerance for each variable in the model
       *
       * Every internal variable receives @p rel_tol as its absolute
       * tolerance.
       *
       * @param[in] rel_tol Solver relative tolerance.
       * @return Zero on success.
       */
      template <typename scalar_type, typename index_type, size_t order>
      int Ieeest<scalar_type, index_type, order>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      /**
       * @brief Residuals of system equations
       *
       * Reads the input signal and evaluates the internal residual.
       *
       * @return Zero on success.
       */
      template <typename scalar_type, typename index_type, size_t order>
      int Ieeest<scalar_type, index_type, order>::evaluateResidual()
      {
        const auto U = static_cast<size_t>(IeeestExternalVariables::U);

        auto* ws = ws_.getData();

        if (auto input = ports_.in.template port<IeeestSignalInputs::input>())
        {
          ws[U]          = input.readSignal();
          ws_indices_[U] = input.signalVariableIndex();
        }

        const auto* y  = y_.getData();
        const auto* yp = yp_.getData();
        auto*       f  = f_.getData();

        evaluateInternalResidual(y, yp, nullptr, ws, f);
        f_.setDataUpdated();
        return 0;
      }

      /**
       * @brief Access the monitor
       *
       * @return Monitor for this model, or nullptr when the model was
       *         constructed without data.
       */
      template <typename scalar_type, typename index_type, size_t order>
      const Model::VariableMonitorBase* Ieeest<scalar_type, index_type, order>::getMonitor() const
      {
        return monitor_.get();
      }

      /**
       * @brief Evaluate the IEEEST internal residual
       *
       * The notch is realized as its two denominator factors in series, the
       * first carrying the numerator. The notch and lead-lag outputs and the
       * unlimited signal are evaluated inline rather than solved for. The
       * notch order is a template parameter, so each specialization is free
       * of branches and loops for sparse automatic differentiation.
       *
       * @param[in] y Internal variables in IeeestInternalVariables order.
       * @param[in] yp Internal derivatives in the same enum order.
       * @param[in] wb Unused; IEEEST has no bus interface.
       * @param[in] ws Signal values in IeeestExternalVariables order.
       * @param[out] f Residuals in IeeestInternalVariables order.
       * @return Zero on success.
       */
      template <typename scalar_type, typename index_type, size_t order>
      __attribute__((always_inline)) inline int
      Ieeest<scalar_type, index_type, order>::evaluateInternalResidual(
          const ScalarT*                  y,
          const ScalarT*                  yp,
          [[maybe_unused]] const ScalarT* wb,
          const ScalarT*                  ws,
          ScalarT*                        f)
      {
        const auto X5  = static_cast<size_t>(IeeestInternalVariables::X5);
        const auto X6  = static_cast<size_t>(IeeestInternalVariables::X6);
        const auto X7  = static_cast<size_t>(IeeestInternalVariables::X7);
        const auto VSS = static_cast<size_t>(IeeestInternalVariables::VSS);
        const auto X1  = static_cast<size_t>(IeeestInternalVariables::X1);
        const auto X2  = static_cast<size_t>(IeeestInternalVariables::X2);
        const auto X3  = static_cast<size_t>(IeeestInternalVariables::X3);
        const auto X4  = static_cast<size_t>(IeeestInternalVariables::X4);
        const auto U   = static_cast<size_t>(IeeestExternalVariables::U);

        const ScalarT x5  = y[X5];
        const ScalarT x6  = y[X6];
        const ScalarT x7  = y[X7];
        const ScalarT vss = y[VSS];

        const ScalarT x5_dot = yp[X5];
        const ScalarT x6_dot = yp[X6];
        const ScalarT x7_dot = yp[X7];

        const ScalarT u = ws[U];

        // Notch output, built factor by factor
        ScalarT v4 = u;
        if constexpr (order == 1)
        {
          const ScalarT x1      = y[X1];
          const ScalarT x1_dot  = yp[X1];
          const ScalarT x1_rate = (u - x1) / b1_;

          f[X1] = -x1_dot + x1_rate;
          v4    = x1 + A5_ * x1_rate;
        }
        if constexpr (order >= 2)
        {
          const ScalarT x1      = y[X1];
          const ScalarT x2      = y[X2];
          const ScalarT x1_dot  = yp[X1];
          const ScalarT x2_dot  = yp[X2];
          const ScalarT x2_rate = (u - x1 - b1_ * x2) / b2_;

          f[X1] = -x1_dot + x2;
          f[X2] = -x2_dot + x2_rate;
          v4    = x1 + A5_ * x2 + A6_ * x2_rate;
        }
        if constexpr (order == 3)
        {
          const ScalarT x3      = y[X3];
          const ScalarT x3_dot  = yp[X3];
          const ScalarT x3_rate = (v4 - x3) / c1_;

          f[X3] = -x3_dot + x3_rate;
          v4    = x3;
        }
        if constexpr (order == 4)
        {
          const ScalarT x3      = y[X3];
          const ScalarT x4      = y[X4];
          const ScalarT x3_dot  = yp[X3];
          const ScalarT x4_dot  = yp[X4];
          const ScalarT x4_rate = (v4 - x3 - c1_ * x4) / c2_;

          f[X3] = -x3_dot + x4;
          f[X4] = -x4_dot + x4_rate;
          v4    = x3;
        }

        // Lead-lags, washout, and output limiter
        const ScalarT x5_rate = (v4 - x5) / T2_;
        const ScalarT v5      = x5 + T1_ * x5_rate;
        const ScalarT x6_rate = (v5 - x6) / T4_;
        const ScalarT v6      = x6 + T3_ * x6_rate;
        const ScalarT x7_rate = (v6 - x7) / T6_;
        const ScalarT v7      = Ks_ * T5_ * x7_rate;

        f[X5]  = -x5_dot + x5_rate;
        f[X6]  = -x6_dot + x6_rate;
        f[X7]  = -x7_dot + x7_rate;
        f[VSS] = -vss + Math::clamp(v7, Lsmin_, Lsmax_);

        return 0;
      }

      //
      //  Private methods
      //

      /**
       * @brief Read the parameters out of the model data
       *
       * No parameter is required; every parameter keeps its documented
       * default when omitted. A non-numeric or non-finite value, or a
       * negative T2, T4, or T6, throws std::invalid_argument. Integer JSON
       * values are accepted for real parameters.
       *
       * @param[in] data Parameters and monitored-variable selections.
       */
      template <typename scalar_type, typename index_type, size_t order>
      void Ieeest<scalar_type, index_type, order>::initializeParameters(const ModelDataT& data)
      {
        using Params = typename ModelDataT::Parameters;

        Model::ParameterReader reader(data, "Ieeest");

        A1_     = reader.readReal(Params::A1).value_or(A1_);
        A2_     = reader.readReal(Params::A2).value_or(A2_);
        A3_     = reader.readReal(Params::A3).value_or(A3_);
        A4_     = reader.readReal(Params::A4).value_or(A4_);
        A5_     = reader.readReal(Params::A5).value_or(A5_);
        A6_     = reader.readReal(Params::A6).value_or(A6_);
        T1_     = reader.readReal(Params::T1).value_or(T1_);
        T2_     = reader.readNonNegative(Params::T2).value_or(T2_);
        T3_     = reader.readReal(Params::T3).value_or(T3_);
        T4_     = reader.readNonNegative(Params::T4).value_or(T4_);
        T5_     = reader.readReal(Params::T5).value_or(T5_);
        T6_     = reader.readNonNegative(Params::T6).value_or(T6_);
        Ks_     = reader.readReal(Params::Ks).value_or(Ks_);
        Lsmin_  = reader.readReal(Params::Lsmin).value_or(Lsmin_);
        Lsmax_  = reader.readReal(Params::Lsmax).value_or(Lsmax_);
        Vcl_    = reader.readReal(Params::Vcl).value_or(Vcl_);
        Vcu_    = reader.readReal(Params::Vcu).value_or(Vcu_);
        Tdelay_ = reader.readReal(Params::Tdelay).value_or(Tdelay_);

        // Model data uses an exact zero to disable these features.
        if (Vcl_ != ZERO<RealT> || Vcu_ != ZERO<RealT> || Tdelay_ != ZERO<RealT>)
        {
          Log::warning() << "Ieeest: Vcl, Vcu, and Tdelay are not modeled and are ignored\n";
        }

        setDerivedParameters();
      }

      /**
       * @brief Bind the monitorable variables to their internal states
       */
      template <typename scalar_type, typename index_type, size_t order>
      void Ieeest<scalar_type, index_type, order>::initializeMonitor()
      {
        using Variable = typename ModelDataT::MonitorableVariables;

        monitor_->set(Variable::vss, [this]
                      { return y_.getData()[static_cast<size_t>(IeeestInternalVariables::VSS)]; });
      }

      /**
       * @brief Log the time-constant floor warning
       *
       * @note Guarded by a static std::atomic<bool> flag so the warning is
       *       printed at most once per model type.
       */
      template <typename scalar_type, typename index_type, size_t order>
      void Ieeest<scalar_type, index_type, order>::logTimeConstantWarning()
      {
        Log::warning() << "Ieeest: T2, T4, and T6 below "
                       << TIME_CONSTANT_MINIMUM
                       << " s are raised to that floor to keep the stabilizer lags well posed\n";
      }

      /**
       * @brief Resolve the parameter-derived constants
       *
       * Raises the lead-lag and washout lags to the well-posedness floor and
       * assigns the notch denominator factors in the order they are realized.
       */
      template <typename scalar_type, typename index_type, size_t order>
      void Ieeest<scalar_type, index_type, order>::setDerivedParameters()
      {
        if (T2_ < TIME_CONSTANT_MINIMUM || T4_ < TIME_CONSTANT_MINIMUM
            || T6_ < TIME_CONSTANT_MINIMUM)
        {
          static std::atomic<bool> time_constant_warning_flag_{false};
          if (!time_constant_warning_flag_.exchange(true))
          {
            logTimeConstantWarning();
          }
        }

        T2_ = std::max(T2_, TIME_CONSTANT_MINIMUM);
        T4_ = std::max(T4_, TIME_CONSTANT_MINIMUM);
        T6_ = std::max(T6_, TIME_CONSTANT_MINIMUM);

        // The factors commute. Up to order two they combine into the first
        // factor; above that the second-order factor goes first to carry the numerator.
        if constexpr (order <= 2)
        {
          b1_ = A1_ + A3_;
          b2_ = A2_ + A4_ + A1_ * A3_;
        }
        else if (factorDegree(A1_, A2_) == 2)
        {
          b1_ = A1_;
          b2_ = A2_;
          c1_ = A3_;
          c2_ = A4_;
        }
        else
        {
          b1_ = A3_;
          b2_ = A4_;
          c1_ = A1_;
          c2_ = A2_;
        }
      }

    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit
