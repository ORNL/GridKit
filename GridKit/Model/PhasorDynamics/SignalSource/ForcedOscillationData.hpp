/**
 * @file ForcedOscillationData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for the forced-oscillation signal source.
 */

#pragma once

#include <cstddef>

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Parameter keys for `ForcedOscillation`; all are optional with documented defaults.
    enum class ForcedOscillationParameters
    {
      A,    ///< \f$A\f$ Oscillation amplitude in output-signal units
      f,    ///< \f$f\f$ Initial oscillation frequency [Hz]
      Kf,   ///< \f$K_f\f$ Linear frequency ramp [Hz/s]
      Phi,  ///< \f$\Phi\f$ Phase offset [rad]
      Ton,  ///< \f$T_{\mathrm{on}}\f$ Activation time [s]
      Toff, ///< \f$T_{\mathrm{off}}\f$ Deactivation time [s]
      Tr,   ///< \f$T_r\f$ Raised-cosine rise time [s]
      Tf,   ///< \f$T_f\f$ Raised-cosine fall time [s]
      Kd    ///< \f$K_d\f$ Exponential decay rate [1/s]
    };

    /// Bus ports for `ForcedOscillation`.
    enum class ForcedOscillationBuses : size_t
    {
      SIZE ///< Number of bus ports; `ForcedOscillation` has no bus attachment
    };

    /// Signal-input ports for `ForcedOscillation`.
    enum class ForcedOscillationSignalInputs : size_t
    {
      SIZE ///< Number of signal-input ports; `ForcedOscillation` has no input
    };

    /// Signal-output ports for `ForcedOscillation`.
    enum class ForcedOscillationSignalOutputs : size_t
    {
      output, ///< \f$s_{\mathrm{FO}}\f$ Required known waveform output
      SIZE    ///< Number of signal-output ports
    };

    /// Variables available through the monitor interface.
    enum class ForcedOscillationMonitorableVariables
    {
      output,   ///< \f$s_{\mathrm{FO}}\f$ Published waveform in output-signal units
      envelope, ///< \f$e\f$ Raised-cosine activation envelope [-]
      active    ///< \f$a\f$ Active-window indicator [-]
    };

    /**
     * @brief Model data for forced-oscillation parameters and signal ports.
     *
     * @tparam real_type Real parameter value type.
     * @tparam index_type Integer index type.
     *
     * @see ForcedOscillation
     */
    template <typename real_type, typename index_type>
    struct ForcedOscillationData : public ComponentData<real_type,
                                                        index_type,
                                                        ForcedOscillationParameters,
                                                        ForcedOscillationBuses,
                                                        ForcedOscillationSignalInputs,
                                                        ForcedOscillationSignalOutputs,
                                                        ForcedOscillationMonitorableVariables>
    {
      ForcedOscillationData() = default;

      using Parameters           = ForcedOscillationParameters;
      using Buses                = ForcedOscillationBuses;
      using SignalInputs         = ForcedOscillationSignalInputs;
      using SignalOutputs        = ForcedOscillationSignalOutputs;
      using MonitorableVariables = ForcedOscillationMonitorableVariables;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
