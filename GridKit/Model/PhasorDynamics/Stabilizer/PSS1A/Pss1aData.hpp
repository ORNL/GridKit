/**
 * @file Pss1aData.hpp
 * @brief Modeling data for the PSS1A stabilizer model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Stabilizer
    {
      /// Parameter keys for the PSS1A stabilizer model.
      enum class Pss1aParameters
      {
        Ks,     ///< \f$K_{\mathrm{S}}\f$ Stabilizer gain
        A1,     ///< \f$A_1\f$ Second-order filter denominator constant
        A2,     ///< \f$A_2\f$ Second-order filter denominator constant
        T1,     ///< \f$T_1\f$ First lead-lag numerator time constant
        T2,     ///< \f$T_2\f$ First lead-lag denominator time constant
        T3,     ///< \f$T_3\f$ Second lead-lag numerator time constant
        T4,     ///< \f$T_4\f$ Second lead-lag denominator time constant
        T5,     ///< \f$T_5\f$ Washout numerator time constant
        T6,     ///< \f$T_6\f$ Input transducer time constant
        Vstmax, ///< \f$V_{\mathrm{STmax}}\f$ Maximum stabilizer output
        Vstmin, ///< \f$V_{\mathrm{STmin}}\f$ Minimum stabilizer output
        Vcu,    ///< \f$V_{\mathrm{cu}}\f$ Upper input cutout threshold
        Vcl     ///< \f$V_{\mathrm{cl}}\f$ Lower input cutout threshold
      };

      /// Buses for the PSS1A stabilizer model.
      enum class Pss1aBuses : size_t
      {
        SIZE
      };

      /// Signal inputs for the PSS1A stabilizer model.
      enum class Pss1aSignalInputs : size_t
      {
        input, ///< Required selected stabilizer input signal ID
        vct,   ///< Optional input-cutout signal ID
        SIZE
      };

      /// Signal outputs for the PSS1A stabilizer model.
      enum class Pss1aSignalOutputs : size_t
      {
        vst, ///< Required stabilizer output signal ID
        SIZE
      };

      /// Variables available through the monitor interface.
      enum class Pss1aMonitorableVariables
      {
        vst ///< Limited stabilizer output
      };

      /**
       * @brief Model data for PSS1A: parameters, signal connections, and
       *        monitored variables.
       *
       * Input selection is represented by explicit signal wiring rather than
       * an integer input-selector parameter.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Pss1aData : public ComponentData<real_type,
                                              index_type,
                                              Pss1aParameters,
                                              Pss1aBuses,
                                              Pss1aSignalInputs,
                                              Pss1aSignalOutputs,
                                              Pss1aMonitorableVariables>
      {
        Pss1aData() = default;

        using Parameters           = Pss1aParameters;
        using Buses                = Pss1aBuses;
        using SignalInputs         = Pss1aSignalInputs;
        using SignalOutputs        = Pss1aSignalOutputs;
        using MonitorableVariables = Pss1aMonitorableVariables;
      };
    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit
