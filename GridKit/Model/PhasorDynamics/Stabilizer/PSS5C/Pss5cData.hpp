/**
 * @file Pss5cData.hpp
 * @brief Modeling data for the PSS5C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Stabilizer
    {
      /// Parameter keys for the PSS5C model.
      enum class Pss5cParameters
      {
        KVL,    ///< \f$K_{\mathrm{VL}}\f$ Very low band gain
        FVL,    ///< \f$F_{\mathrm{VL}}\f$ Very low band central frequency
        VVLmax, ///< \f$V_{\mathrm{VLmax}}\f$ Very low band upper limit
        VVLmin, ///< \f$V_{\mathrm{VLmin}}\f$ Very low band lower limit
        KL,     ///< \f$K_{\mathrm{L}}\f$ Low band gain
        FL,     ///< \f$F_{\mathrm{L}}\f$ Low band central frequency
        VLmax,  ///< \f$V_{\mathrm{Lmax}}\f$ Low band upper limit
        VLmin,  ///< \f$V_{\mathrm{Lmin}}\f$ Low band lower limit
        KI,     ///< \f$K_{\mathrm{I}}\f$ Intermediate band gain
        FI,     ///< \f$F_{\mathrm{I}}\f$ Intermediate band central frequency
        VImax,  ///< \f$V_{\mathrm{Imax}}\f$ Intermediate band upper limit
        VImin,  ///< \f$V_{\mathrm{Imin}}\f$ Intermediate band lower limit
        KH,     ///< \f$K_{\mathrm{H}}\f$ High band gain
        FH,     ///< \f$F_{\mathrm{H}}\f$ High band central frequency
        VHmax,  ///< \f$V_{\mathrm{Hmax}}\f$ High band upper limit
        VHmin,  ///< \f$V_{\mathrm{Hmin}}\f$ High band lower limit
        k1,     ///< \f$k_{1}\f$
        k2,     ///< \f$k_{2}\f$
        k3,     ///< \f$k_{3}\f$
        VSTMAX, ///< \f$V_{\mathrm{STMAX}}\f$ Maximum PSS output
        VSTMIN  ///< \f$V_{\mathrm{STMIN}}\f$ Minimum PSS output
      };

      /// Bus keys for the PSS5C model; deferred until port integration.
      enum class Pss5cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the PSS5C model; deferred until port integration.
      enum class Pss5cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the PSS5C model; deferred until port integration.
      enum class Pss5cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the PSS5C model; deferred until implementation.
      enum class Pss5cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for PSS5C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Pss5cData : public ComponentData<real_type,
                                              index_type,
                                              Pss5cParameters,
                                              Pss5cBuses,
                                              Pss5cSignalInputs,
                                              Pss5cSignalOutputs,
                                              Pss5cMonitorableVariables>
      {
        Pss5cData() = default;

        using Parameters           = Pss5cParameters;
        using Buses                = Pss5cBuses;
        using SignalInputs         = Pss5cSignalInputs;
        using SignalOutputs        = Pss5cSignalOutputs;
        using MonitorableVariables = Pss5cMonitorableVariables;
      };
    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit
