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
      /// Parameter keys for the ESDC1A exciter model. Every parameter is
      /// optional and retains its documented default when omitted.
      enum class Esdc1aParameters
      {
        Tr,     ///< \f$T_R\f$ Voltage transducer time constant [sec]
        Ka,     ///< \f$K_A\f$ Voltage-regulator gain [p.u.]
        Ta,     ///< \f$T_A\f$ Voltage-regulator time constant [sec]
        Tb,     ///< \f$T_B\f$ Input lead-lag denominator time constant [sec]
        Tc,     ///< \f$T_C\f$ Input lead-lag numerator time constant [sec]
        Vrmax,  ///< \f$V_R^{\max}\f$ Maximum voltage-regulator output [p.u.]
        Vrmin,  ///< \f$V_R^{\min}\f$ Minimum voltage-regulator output [p.u.]
        Ke,     ///< \f$K_E\f$ Exciter constant [p.u.]
        Te,     ///< \f$T_E\f$ Exciter time constant [sec]
        Kf,     ///< \f$K_F\f$ Stabilizing feedback gain [p.u.]
        Tf1,    ///< \f$T_{F1}\f$ Stabilizing feedback time constant [sec]
        Spdmlt, ///< \f$s_{\mathrm{spd}}\f$ Field-voltage speed-multiplier flag [boolean]
        E1,     ///< \f$E_1\f$ First saturation voltage point [p.u.]
        Se1,    ///< \f$S_E(E_1)\f$ Saturation coefficient at \f$E_1\f$ [p.u.]
        E2,     ///< \f$E_2\f$ Second saturation voltage point [p.u.]
        Se2,    ///< \f$S_E(E_2)\f$ Saturation coefficient at \f$E_2\f$ [p.u.]
        UEL,    ///< \f$I_{\mathrm{UEL}}\f$ UEL input-routing selector [integer]
        exclim  ///< \f$s_{\mathrm{lim}}\f$ Exciter field-voltage-state lower-limit flag [boolean]
      };

      /// Buses for the ESDC1A exciter model.
      enum class Esdc1aBuses : size_t
      {
        bus, ///< \f$V_{\mathrm{r}},V_{\mathrm{i}}\f$ Required Known terminal-bus voltage [p.u.]
        SIZE ///< Number of ESDC1A bus ports
      };

      /// Signal inputs for the ESDC1A exciter model.
      enum class Esdc1aSignalInputs : size_t
      {
        speed, ///< \f$\omega\f$ Known machine speed-deviation input [p.u.]; required when \f$s_{\mathrm{spd}}=1\f$
        vref,  ///< \f$V_{\mathrm{ref}}\f$ Optional Unknown voltage-reference input [p.u.]
        vs,    ///< \f$V_S\f$ Optional Known stabilizer input [p.u.]
        vuel,  ///< \f$V_{\mathrm{UEL}}\f$ Optional Known UEL input [p.u.]
        SIZE   ///< Number of ESDC1A input-signal ports
      };

      /// Signal outputs for the ESDC1A exciter model.
      enum class Esdc1aSignalOutputs : size_t
      {
        efd, ///< \f$E_{\mathrm{fd}}\f$ Required Known field-voltage output [p.u.]
        SIZE ///< Number of ESDC1A output-signal ports
      };

      /// Variables available through the monitor interface.
      enum class Esdc1aMonitorableVariables
      {
        efd, ///< \f$E_{\mathrm{fd}}\f$ Field-voltage output [p.u.]
        vc,  ///< \f$V_C\f$ Filtered terminal-voltage magnitude [p.u.]
        vr,  ///< \f$V_R\f$ Voltage-regulator output [p.u.]
        vf,  ///< \f$V_F\f$ Stabilizing feedback state [p.u.]
        se,  ///< \f$E_{\mathrm{fd}}'S_E(E_{\mathrm{fd}}')\f$ Scaled-quadratic saturation contribution [p.u.]
        vfe  ///< \f$V_{\mathrm{FE}}\f$ Exciter feedback drive [p.u.]
      };

      /**
       * @brief Model data for ESDC1A parameters, terminal bus, signal ports,
       *        and monitored variables.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       *
       * @see Esdc1a
       */
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
