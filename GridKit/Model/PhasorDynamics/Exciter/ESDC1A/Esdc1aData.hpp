/**
 * @file Esdc1aData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for the ESDC1A exciter model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Parameter keys for the ESDC1A exciter model.
      enum class Esdc1aParameters
      {
        Tr,     ///< \f$T_R\f$ Voltage transducer time constant
        Ka,     ///< \f$K_A\f$ Voltage-regulator gain
        Ta,     ///< \f$T_A\f$ Voltage-regulator time constant
        Tb,     ///< \f$T_B\f$ Input lead-lag denominator time constant
        Tc,     ///< \f$T_C\f$ Input lead-lag numerator time constant
        Vrmax,  ///< \f$V_R^{\max}\f$ Maximum voltage-regulator output
        Vrmin,  ///< \f$V_R^{\min}\f$ Minimum voltage-regulator output
        Ke,     ///< \f$K_E\f$ Exciter constant
        Te,     ///< \f$T_E\f$ Exciter time constant
        Kf,     ///< \f$K_F\f$ Stabilizing feedback gain
        Tf1,    ///< \f$T_{F1}\f$ Stabilizing feedback time constant
        Spdmlt, ///< \f$s_{\mathrm{spd}}\f$ Field-voltage speed-multiplier flag
        E1,     ///< \f$E_1\f$ First saturation voltage point
        Se1,    ///< \f$S_E(E_1)\f$ Saturation coefficient at \f$E_1\f$
        E2,     ///< \f$E_2\f$ Second saturation voltage point
        Se2,    ///< \f$S_E(E_2)\f$ Saturation coefficient at \f$E_2\f$
        UEL,    ///< \f$I_{\mathrm{UEL}}\f$ UEL input-routing selector
        exclim  ///< \f$s_{\mathrm{lim}}\f$ Exciter feedback lower-limit flag
      };

      /// Buses for the ESDC1A exciter model.
      enum class Esdc1aBuses : size_t
      {
        bus, ///< Terminal bus ID for \f$V_{\mathrm{r}}\f$ and \f$V_{\mathrm{i}}\f$
        SIZE
      };

      /// Signal inputs for the ESDC1A exciter model.
      enum class Esdc1aSignalInputs : size_t
      {
        speed, ///< \f$\omega\f$ Machine speed-deviation signal ID
        vref,  ///< \f$V_{\mathrm{ref}}\f$ Optional voltage-reference signal ID
        vs,    ///< \f$V_S\f$ Optional stabilizer input signal ID
        vuel,  ///< \f$V_{\mathrm{UEL}}\f$ Optional UEL input signal ID
        SIZE
      };

      /// Signal outputs for the ESDC1A exciter model.
      enum class Esdc1aSignalOutputs : size_t
      {
        efd, ///< \f$E_{\mathrm{fd}}\f$ Required field-voltage output signal ID
        SIZE
      };

      /// Variables available through the monitor interface.
      enum class Esdc1aMonitorableVariables
      {
        efd, ///< \f$E_{\mathrm{fd}}\f$ Field-voltage output
        vc,  ///< \f$V_C\f$ Filtered terminal-voltage magnitude
        vr,  ///< \f$V_R\f$ Voltage-regulator output
        vf,  ///< \f$V_F\f$ Stabilizing feedback state
        se,  ///< \f$S_E\f$ Exciter saturation coefficient
        vfe  ///< \f$V_{\mathrm{FE}}\f$ Exciter feedback drive
      };

      template <typename real_type, typename index_type>
      struct Esdc1aData : public ComponentData<real_type,
                                               index_type,
                                               Esdc1aParameters,
                                               Esdc1aBuses,
                                               Esdc1aSignalInputs,
                                               Esdc1aSignalOutputs,
                                               Esdc1aMonitorableVariables>
      {
        Esdc1aData() = default;

        using Parameters           = Esdc1aParameters;
        using Buses                = Esdc1aBuses;
        using SignalInputs         = Esdc1aSignalInputs;
        using SignalOutputs        = Esdc1aSignalOutputs;
        using MonitorableVariables = Esdc1aMonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
