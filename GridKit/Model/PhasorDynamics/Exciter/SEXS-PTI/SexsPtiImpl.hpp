#pragma once

/**
 * @file SexsPtiImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the SEXS-PTI exciter model.
 */

#include <algorithm>
#include <atomic>
#include <format>
#include <iostream>

#include <GridKit/Model/ConfigurationChecks.hpp>
#include <GridKit/Model/ParameterReader.hpp>
#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/SEXS-PTI/SexsPti.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/SEXS-PTI/SexsPtiData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeSet.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      using Log = ::GridKit::Utilities::Logger;

      template <typename scalar_type, typename index_type>
      SexsPti<scalar_type, index_type>::SexsPti(BusT* bus)
        : bus_(bus)
      {
        size_ = 3;
        setDerivedParameters();
      }

      template <typename scalar_type, typename index_type>
      SexsPti<scalar_type, index_type>::SexsPti(BusT* bus, const ModelDataT& data)
        : bus_(bus),
          monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        initializeMonitor();
        size_ = 3;
      }

      template <typename scalar_type, typename index_type>
      SexsPti<scalar_type, index_type>::~SexsPti()
      {
      }

      template <typename scalar_type, typename index_type>
      int SexsPti<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int SexsPti<scalar_type, index_type>::allocate()
      {
        if (!allocated_)
        {
          this->allocateVectors(size_);
        }
        auto size = static_cast<size_t>(size_);

        tag_.resize(size);

        variable_indices_.resize(size);
        residual_indices_.resize(size);

        // Default variable and residual index mapping to local index
        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        wb_.resize(2);

        const auto signal_size = Utilities::enum_size<SexsPtiExternalVariables>();
        ws_.resize(static_cast<IdxT>(signal_size));
        ws_.setToZero();
        ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

        if (auto efd_port = ports_.out.template port<SexsPtiSignalOutputs::efd>())
        {
          efd_port.link(&y_.getData()[1], &(this->getVariableIndex(1)));
        }

        allocated_ = true;
        return 0;
      }

      /**
       * @brief Validate the SEXS-PTI configuration
       *
       * Checks the terminal bus, the regulator gain initialization divides
       * by, the field-voltage limits, the required field-voltage output, and
       * attached external signals.
       *
       * @return The configuration checks; passed() when valid.
       */
      template <typename scalar_type, typename index_type>
      Model::ConfigurationChecks SexsPti<scalar_type, index_type>::verify() const
      {
        Model::ConfigurationChecks checks("SexsPti");

        checks.check(bus_ != nullptr, "bus pointer is null");
        checks.check(K_ > ZERO<RealT>, "K must be positive");
        checks.check(Efdmin_ < Efdmax_, "Efdmin must be less than Efdmax");
        checks.check(ports_.out.template port<SexsPtiSignalOutputs::efd>().connected(),
                     "required efd output signal is not assigned");

        // An attached port must resolve to readable signal storage. The
        // enumerator is a template argument, so each port names itself once.
        auto check_attached_signal =
            [&]<SexsPtiSignalInputs input>(const char* name)
        {
          if (ports_.in.template port<input>().connected()
              && !ports_.in.template port<input>().linked())
          {
            checks.fail(std::format("{} signal attached with no linked source", name));
          }
        };

        check_attached_signal.template operator()<SexsPtiSignalInputs::vref>("vref");
        check_attached_signal.template operator()<SexsPtiSignalInputs::vs>("vs");
        check_attached_signal.template operator()<SexsPtiSignalInputs::vuel>("vuel");
        check_attached_signal.template operator()<SexsPtiSignalInputs::voel>("voel");

        return checks;
      }

      template <typename scalar_type, typename index_type>
      int SexsPti<scalar_type, index_type>::initialize()
      {
        if (!this->verifyAndLog())
        {
          Log::error() << "SexsPti: cannot initialize with invalid configuration\n";
          return 1;
        }

        auto* y  = y_.getData();
        auto* yp = yp_.getData();

        // The assigned efd node aliases this entry after allocate(), so it
        // holds the field voltage seeded by the machine.
        const ScalarT efd0 = y[1];

        // Setpoint members provide the defaults for unattached signals.
        auto read_signal = [&]<SexsPtiSignalInputs input>(const ScalarT& default_value) -> ScalarT
        {
          if (auto port = ports_.in.template port<input>())
          {
            return port.readSignal();
          }
          return default_value;
        };

        const ScalarT vs   = read_signal.template operator()<SexsPtiSignalInputs::vs>(vs_set_);
        const ScalarT vuel = read_signal.template operator()<SexsPtiSignalInputs::vuel>(vuel_set_);
        const ScalarT voel = read_signal.template operator()<SexsPtiSignalInputs::voel>(voel_set_);

        uel_on_ = ZERO<RealT>;
        if (ports_.in.template port<SexsPtiSignalInputs::vuel>())
        {
          uel_on_ = ONE<RealT>;
        }

        oel_on_ = ZERO<RealT>;
        if (ports_.in.template port<SexsPtiSignalInputs::voel>())
        {
          oel_on_ = ONE<RealT>;
        }

        ScalarT vreal = bus_->Vr();
        ScalarT vimag = bus_->Vi();
        ScalarT Ec    = std::sqrt(vreal * vreal + vimag * vimag);
        ScalarT vtr   = efd0 / K_;
        ScalarT vr    = (Ta_ - Tb_) * vtr;
        ScalarT vref  = Ec + vtr - vs - uel_on_ * vuel - oel_on_ * voel;

        y[0] = vr;
        y[1] = efd0;
        y[2] = vtr;

        for (IdxT i = 0; i < size_; ++i)
        {
          yp[static_cast<size_t>(i)] = 0.0;
        }

        vref_set_ = static_cast<RealT>(vref);
        vs_set_   = static_cast<RealT>(vs);
        vuel_set_ = static_cast<RealT>(vuel);
        voel_set_ = static_cast<RealT>(voel);

        if (auto vref_port = ports_.in.template port<SexsPtiSignalInputs::vref>())
        {
          vref_port.writeValue(vref_set_);
        }

        y_.setDataUpdated();
        yp_.setDataUpdated();

        // For DependencyTracking::Variable, set variable numbers
        if constexpr (std::is_same_v<scalar_type, DependencyTracking::Variable>)
        {
          this->initializeDependencyTrackingVariableNumbers();
        }

        return 0;
      }

      template <typename scalar_type, typename index_type>
      int SexsPti<scalar_type, index_type>::tagDifferentiable()
      {
        tag_[0] = true;
        tag_[1] = true;
        tag_[2] = false;

        return 0;
      }

      /**
       * @brief Compute the absolute tolerance for each variable in the model
       *
       * @param rel_tol The relative tolerance which can be used to pick the
       *        absolute tolerance.
       * @tparam scalar_type Scalar data type
       * @tparam index_type Index data type
       * @return int 0 if successful, non-zero otherwise.
       *
       * This represents a "noise" level close to zero for which pure relative
       * error cannot be used.
       */
      template <typename scalar_type, typename index_type>
      int SexsPti<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline int SexsPti<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT* y,
          const ScalarT* yp,
          const ScalarT* wb,
          const ScalarT* ws,
          ScalarT*       f)
      {
        const auto VREF = static_cast<size_t>(SexsPtiExternalVariables::VREF);
        const auto VS   = static_cast<size_t>(SexsPtiExternalVariables::VS);
        const auto VUEL = static_cast<size_t>(SexsPtiExternalVariables::VUEL);
        const auto VOEL = static_cast<size_t>(SexsPtiExternalVariables::VOEL);

        ScalarT vr      = y[0];
        ScalarT efd     = y[1];
        ScalarT vtr     = y[2];
        ScalarT vr_dot  = yp[0];
        ScalarT efd_dot = yp[1];

        ScalarT Ec   = std::sqrt(wb[0] * wb[0] + wb[1] * wb[1]);
        ScalarT vref = ws[VREF];
        ScalarT vs   = ws[VS];
        ScalarT vuel = ws[VUEL];
        ScalarT voel = ws[VOEL];

        ScalarT func = (-efd + (K_ / Tb_) * (-vr + Ta_ * vtr)) / Te_;

        f[0] = -vr_dot + (-vr + Ta_ * vtr) / Tb_ - vtr;
        f[1] = -efd_dot + Math::antiwindup(efd, func, Efdmin_, Efdmax_);
        f[2] = -vtr - Ec + vref + vs + uel_on_ * vuel + oel_on_ * voel;

        return 0;
      }

      template <typename scalar_type, typename index_type>
      int SexsPti<scalar_type, index_type>::evaluateResidual()
      {
        auto* ws = ws_.getData();

        // Attached signals are read live; unattached ones keep the latched value.
        auto read_signal = [&]<SexsPtiSignalInputs      input,
                               SexsPtiExternalVariables variable>(const ScalarT& latched)
        {
          const auto index   = static_cast<size_t>(variable);
          ws[index]          = latched;
          ws_indices_[index] = INVALID_INDEX<IdxT>;
          if (auto port = ports_.in.template port<input>())
          {
            ws[index]          = port.readSignal();
            ws_indices_[index] = port.signalVariableIndex();
          }
        };

        read_signal.template operator()<SexsPtiSignalInputs::vref,
                                        SexsPtiExternalVariables::VREF>(vref_set_);
        read_signal.template operator()<SexsPtiSignalInputs::vs,
                                        SexsPtiExternalVariables::VS>(vs_set_);
        read_signal.template operator()<SexsPtiSignalInputs::vuel,
                                        SexsPtiExternalVariables::VUEL>(vuel_set_);
        read_signal.template operator()<SexsPtiSignalInputs::voel,
                                        SexsPtiExternalVariables::VOEL>(voel_set_);

        auto* wb = wb_.getData();
        wb[0]    = bus_->Vr();
        wb[1]    = bus_->Vi();

        const auto* y  = y_.getData();
        const auto* yp = yp_.getData();
        auto*       f  = f_.getData();
        evaluateInternalResidual(y, yp, wb, ws, f);

        f_.setDataUpdated();

        return 0;
      }

      /**
       * @brief Read the parameters out of the model data
       *
       * Every parameter is optional and keeps its documented default when
       * omitted. A non-numeric or non-finite value, or a negative floored
       * time constant, throws std::invalid_argument. Integer JSON values are
       * accepted for real parameters.
       *
       * @param[in] data Parameters and monitored-variable selections.
       */
      template <typename scalar_type, typename index_type>
      void SexsPti<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Params = typename ModelDataT::Parameters;

        const Model::ParameterReader reader(data, "SexsPti");

        Ta_     = reader.readReal(Params::Ta).value_or(Ta_);
        Tb_     = reader.readNonNegative(Params::Tb).value_or(Tb_);
        Te_     = reader.readNonNegative(Params::Te).value_or(Te_);
        K_      = reader.readReal(Params::K).value_or(K_);
        Efdmax_ = reader.readReal(Params::Efdmax).value_or(Efdmax_);
        Efdmin_ = reader.readReal(Params::Efdmin).value_or(Efdmin_);
        setDerivedParameters();
      }

      /**
       * @brief Static method to log time constant warnings
       *
       * @note Guarded by a static std::atomic<bool> flag so the warning is
       *       printed at most once per model type.
       */
      template <typename scalar_type, typename index_type>
      void SexsPti<scalar_type, index_type>::logTimeConstantWarning()
      {
        Log::warning() << "SexsPti: Tb and Te below " << TIME_CONSTANT_MINIMUM
                       << " s are raised to that floor\n";
      }

      /**
       * @brief Resolve the parameter-derived constants
       *
       * Raises the lead-lag denominator and exciter lags to the
       * well-posedness floor.
       */
      template <typename scalar_type, typename index_type>
      void SexsPti<scalar_type, index_type>::setDerivedParameters()
      {
        if (Tb_ < TIME_CONSTANT_MINIMUM || Te_ < TIME_CONSTANT_MINIMUM)
        {
          static std::atomic<bool> time_constant_warning_flag_{false};
          if (!time_constant_warning_flag_.exchange(true))
          {
            logTimeConstantWarning();
          }
        }

        Tb_ = std::max(Tb_, TIME_CONSTANT_MINIMUM);
        Te_ = std::max(Te_, TIME_CONSTANT_MINIMUM);
      }

      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* SexsPti<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      template <typename scalar_type, typename index_type>
      void SexsPti<scalar_type, index_type>::initializeMonitor()
      {
        using Variable = typename ModelDataT::MonitorableVariables;
        monitor_->set(Variable::efd, [this]
                      { return y_.getData()[1]; });
      }

    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
