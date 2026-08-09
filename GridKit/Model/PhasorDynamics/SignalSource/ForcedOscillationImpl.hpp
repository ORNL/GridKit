/**
 * @file ForcedOscillationImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the forced-oscillation signal source.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <variant>

#include <GridKit/Model/PhasorDynamics/SignalSource/ForcedOscillation.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Logger used for `ForcedOscillation` diagnostics.
    using Log = ::GridKit::Utilities::Logger;

    /**
     * @brief Construct an unconfigured zero-amplitude source.
     */
    template <typename scalar_type, typename index_type>
    ForcedOscillation<scalar_type, index_type>::ForcedOscillation()
    {
      size_ = 0;
      refreshOutput(static_cast<RealT>(0.0));
    }

    /**
     * @brief Construct a source from model data.
     *
     * @param[in] data Parameters and signal-port selections.
     */
    template <typename scalar_type, typename index_type>
    ForcedOscillation<scalar_type, index_type>::ForcedOscillation(const ModelDataT& data)
      : monitor_(std::make_unique<MonitorT>(data))
    {
      initializeParameters(data);
      initializeMonitor();
      size_ = 0;
      refreshOutput(static_cast<RealT>(0.0));
    }

    /**
     * @brief Destroy the stateless source.
     */
    template <typename scalar_type, typename index_type>
    ForcedOscillation<scalar_type, index_type>::~ForcedOscillation()
    {
    }

    /**
     * @brief Set the component identifier assigned by the system model.
     *
     * @param[in] component_id System-model component identifier.
     */
    template <typename scalar_type, typename index_type>
    int ForcedOscillation<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /**
     * @brief Bind the output signal to source-owned storage.
     *
     * The output deliberately retains `INVALID_INDEX` because it is an
     * exogenous value rather than a system DAE variable.
     */
    template <typename scalar_type, typename index_type>
    int ForcedOscillation<scalar_type, index_type>::allocate()
    {
      static constexpr auto OUTPUT = ForcedOscillationInternalVariables::OUTPUT;

      output_index_ = INVALID_INDEX<IdxT>;
      if (signals_.template isAssigned<OUTPUT>())
      {
        signals_.template getSignalNode<OUTPUT>()->set(&output_, &output_index_);
      }

      allocated_ = true;
      return 0;
    }

    /**
     * @brief Validate parameters and the required output assignment.
     *
     * @return Number of configuration errors; zero when valid.
     */
    template <typename scalar_type, typename index_type>
    int ForcedOscillation<scalar_type, index_type>::verify() const
    {
      static constexpr auto OUTPUT = ForcedOscillationInternalVariables::OUTPUT;

      int ret = static_cast<int>(parameter_error_count_);

      auto check = [&](bool condition, const char* message)
      {
        if (!condition)
        {
          Log::error() << "ForcedOscillation: " << message << '\n';
          ++ret;
        }
      };

      check(std::isfinite(A_) && A_ >= static_cast<RealT>(0.0),
            "A must be finite and non-negative");
      check(std::isfinite(f_) && f_ >= static_cast<RealT>(0.0),
            "f must be finite and non-negative");
      check(std::isfinite(Kf_) && Kf_ >= static_cast<RealT>(0.0),
            "Kf must be finite and non-negative");
      check(std::isfinite(Phi_), "Phi must be finite");
      check(std::isfinite(Ton_), "Ton must be finite");
      check(std::isfinite(Toff_), "Toff must be finite");
      if (std::isfinite(Ton_) && std::isfinite(Toff_))
      {
        check(Toff_ < static_cast<RealT>(0.0) || Toff_ >= Ton_,
              "Toff must be negative or greater than or equal to Ton");
      }
      check(std::isfinite(Tr_) && Tr_ >= static_cast<RealT>(0.0),
            "Tr must be finite and non-negative");
      check(std::isfinite(Tf_) && Tf_ >= static_cast<RealT>(0.0),
            "Tf must be finite and non-negative");
      check(std::isfinite(Kd_) && Kd_ >= static_cast<RealT>(0.0),
            "Kd must be finite and non-negative");
      check(signals_.template isAssigned<OUTPUT>(),
            "required output signal is not assigned");

      return ret;
    }

    /**
     * @brief Initialize the published waveform at the current model time.
     */
    template <typename scalar_type, typename index_type>
    int ForcedOscillation<scalar_type, index_type>::initialize()
    {
      refreshOutput(time_);
      return 0;
    }

    /**
     * @brief Tag differential variables.
     *
     * This source owns no system variables.
     */
    template <typename scalar_type, typename index_type>
    int ForcedOscillation<scalar_type, index_type>::tagDifferentiable()
    {
      return 0;
    }

    /**
     * @brief Set absolute tolerances.
     *
     * This source owns no system variables.
     */
    template <typename scalar_type, typename index_type>
    int ForcedOscillation<scalar_type, index_type>::setAbsoluteTolerance(RealT)
    {
      return 0;
    }

    /**
     * @brief Refresh the source value for the current solver time.
     *
     * @param[in] time Current simulation time [s].
     * @param[in] alpha Solver Jacobian shift.
     */
    template <typename scalar_type, typename index_type>
    void ForcedOscillation<scalar_type, index_type>::updateTime(RealT time, RealT alpha)
    {
      Component<scalar_type, index_type>::updateTime(time, alpha);
      refreshOutput(time);
    }

    /**
     * @brief Evaluate residuals.
     *
     * This source owns no residual rows.
     */
    template <typename scalar_type, typename index_type>
    int ForcedOscillation<scalar_type, index_type>::evaluateResidual()
    {
      return 0;
    }

    /**
     * @brief Construct the empty Jacobian for this stateless source.
     */
    template <typename scalar_type, typename index_type>
    int ForcedOscillation<scalar_type, index_type>::evaluateJacobian()
    {
      return this->constructCoo();
    }

    /**
     * @brief Return the source's signal interface.
     */
    template <typename scalar_type, typename index_type>
    auto ForcedOscillation<scalar_type, index_type>::getSignals()
        -> ComponentSignals<ScalarT,
                            IdxT,
                            ForcedOscillationInternalVariables,
                            ForcedOscillationExternalVariables>&
    {
      return signals_;
    }

    /**
     * @brief Return the configured variable monitor, if any.
     */
    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase*
    ForcedOscillation<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    //
    // Private methods
    //

    /**
     * @brief Load one optional real-valued parameter.
     *
     * Real and integer serialized values are accepted. Any other stored type
     * records a loading error while preserving the documented default.
     */
    template <typename scalar_type, typename index_type>
    void ForcedOscillation<scalar_type, index_type>::loadRealParameter(
        const ModelDataT&           data,
        ForcedOscillationParameters parameter,
        RealT&                      target,
        const char*                 name)
    {
      if (!data.parameters.contains(parameter))
      {
        return;
      }

      const auto& value = data.parameters.at(parameter);
      if (const auto* real_value = std::get_if<RealT>(&value))
      {
        target = *real_value;
      }
      else if (const auto* index_value = std::get_if<IdxT>(&value))
      {
        target = static_cast<RealT>(*index_value);
      }
      else
      {
        Log::error() << "ForcedOscillation: parameter '" << name
                     << "' must be numeric\n";
        ++parameter_error_count_;
      }
    }

    /**
     * @brief Read optional parameters while preserving omitted-key defaults.
     */
    template <typename scalar_type, typename index_type>
    void ForcedOscillation<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Parameters = typename ModelDataT::Parameters;

      parameter_error_count_ = 0;

      loadRealParameter(data, Parameters::A, A_, "A");
      loadRealParameter(data, Parameters::f, f_, "f");
      loadRealParameter(data, Parameters::Kf, Kf_, "Kf");
      loadRealParameter(data, Parameters::Phi, Phi_, "Phi");
      loadRealParameter(data, Parameters::Ton, Ton_, "Ton");
      loadRealParameter(data, Parameters::Toff, Toff_, "Toff");
      loadRealParameter(data, Parameters::Tr, Tr_, "Tr");
      loadRealParameter(data, Parameters::Tf, Tf_, "Tf");
      loadRealParameter(data, Parameters::Kd, Kd_, "Kd");
    }

    /**
     * @brief Bind monitor selections to the published waveform diagnostics.
     */
    template <typename scalar_type, typename index_type>
    void ForcedOscillation<scalar_type, index_type>::initializeMonitor()
    {
      using Variable = typename ModelDataT::MonitorableVariables;

      monitor_->set(Variable::output, [this]
                    { return output_; });
      monitor_->set(Variable::envelope, [this]
                    { return envelope_; });
      monitor_->set(Variable::active, [this]
                    { return active_; });
    }

    /**
     * @brief Evaluate and publish the configured waveform.
     *
     * @param[in] time Simulation time [s].
     */
    template <typename scalar_type, typename index_type>
    void ForcedOscillation<scalar_type, index_type>::refreshOutput(RealT time)
    {
      static constexpr RealT PI = static_cast<RealT>(
          3.141592653589793238462643383279502884L);

      const bool after_start = time >= Ton_;
      const bool before_stop = Toff_ < static_cast<RealT>(0.0) || time < Toff_;
      const bool active      = after_start && before_stop;

      RealT envelope = static_cast<RealT>(0.0);
      if (active)
      {
        envelope = static_cast<RealT>(1.0);

        if (Tr_ > static_cast<RealT>(0.0) && time < Ton_ + Tr_)
        {
          const RealT x  = (time - Ton_) / Tr_;
          envelope      *= static_cast<RealT>(0.5)
                      * (static_cast<RealT>(1.0) - std::cos(PI * x));
        }

        if (Toff_ >= static_cast<RealT>(0.0)
            && Tf_ > static_cast<RealT>(0.0)
            && time > Toff_ - Tf_)
        {
          const RealT x  = (Toff_ - time) / Tf_;
          envelope      *= static_cast<RealT>(0.5)
                      * (static_cast<RealT>(1.0) - std::cos(PI * x));
        }
      }

      const RealT tau   = std::max(time - Ton_, static_cast<RealT>(0.0));
      const RealT phase = Phi_ + static_cast<RealT>(2.0) * PI * (f_ * tau + static_cast<RealT>(0.5) * Kf_ * tau * tau);
      const RealT decay = std::exp(-Kd_ * tau);

      envelope_ = envelope;
      active_   = active ? static_cast<RealT>(1.0) : static_cast<RealT>(0.0);
      output_   = A_ * envelope * decay * std::sin(phase);
    }
  } // namespace PhasorDynamics
} // namespace GridKit
