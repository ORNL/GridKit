/**
 * @file Uel1Data.hpp
 * @brief Modeling data for the UEL1 model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Limiter
    {
      /// Parameter keys for the UEL1 model.
      enum class Uel1Parameters
      {
        KUC,     ///< \f$K_{\mathrm{UC}}\f$ UEL center setting
        KUR,     ///< \f$K_{\mathrm{UR}}\f$ UEL radius setting
        VURmax,  ///< \f$V_{\mathrm{URmax}}\f$
        VUCmax,  ///< \f$V_{\mathrm{UCmax}}\f$
        KUF,     ///< \f$K_{\mathrm{UF}}\f$ UEL excitation system stabilizer gain
        KUI,     ///< \f$K_{\mathrm{UI}}\f$ UEL integral gain
        KUL,     ///< \f$K_{\mathrm{UL}}\f$ UEL proportional gain
        VUImax,  ///< \f$V_{\mathrm{UImax}}\f$ UEL PI control maximum output
        VUImin,  ///< \f$V_{\mathrm{UImin}}\f$ UEL PI control minimum output
        TU1,     ///< \f$T_{\mathrm{U1}}\f$ UEL numerator (lead) time constant (first block)
        TU2,     ///< \f$T_{\mathrm{U2}}\f$ UEL denominator (lag) time constant (first block)
        TU3,     ///< \f$T_{\mathrm{U3}}\f$ UEL numerator (lead) time constant (second block)
        TU4,     ///< \f$T_{\mathrm{U4}}\f$ UEL denominator (lag) time constant (second block)
        VUELmax, ///< \f$V_{\mathrm{UELmax}}\f$ UEL maximum output
        VUELmin  ///< \f$V_{\mathrm{UELmin}}\f$ UEL minimum output
      };

      /// Bus keys for the UEL1 model; deferred until port integration.
      enum class Uel1Buses : size_t
      {
        SIZE
      };

      /// Signal input keys for the UEL1 model; deferred until port integration.
      enum class Uel1SignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the UEL1 model; deferred until port integration.
      enum class Uel1SignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the UEL1 model; deferred until implementation.
      enum class Uel1MonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for UEL1.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Uel1Data : public ComponentData<real_type,
                                             index_type,
                                             Uel1Parameters,
                                             Uel1Buses,
                                             Uel1SignalInputs,
                                             Uel1SignalOutputs,
                                             Uel1MonitorableVariables>
      {
        Uel1Data() = default;

        using Parameters           = Uel1Parameters;
        using Buses                = Uel1Buses;
        using SignalInputs         = Uel1SignalInputs;
        using SignalOutputs        = Uel1SignalOutputs;
        using MonitorableVariables = Uel1MonitorableVariables;
      };
    } // namespace Limiter
  } // namespace PhasorDynamics
} // namespace GridKit
