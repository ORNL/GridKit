/**
 * @file Pss3cData.hpp
 * @brief Modeling data for the PSS3C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Stabilizer
    {
      /// Parameter keys for the PSS3C model.
      enum class Pss3cParameters
      {
        KS1,    ///< \f$K_{\mathrm{S1}}\f$ Power system stabilizer gain for input channel 1
        KS2,    ///< \f$K_{\mathrm{S2}}\f$ Power system stabilizer gain for input channel 2
        T1,     ///< \f$T_{1}\f$ PSS transducer time constant for input channel 1
        T2,     ///< \f$T_{2}\f$ PSS transducer time constant for input channel 2
        Tw1,    ///< \f$T_{\mathrm{w1}}\f$ washout time constant (input channel 1)
        Tw2,    ///< \f$T_{\mathrm{w2}}\f$ washout time constant (input channel 2)
        Tw3,    ///< \f$T_{\mathrm{w3}}\f$ washout time constant (combined channels)
        A1,     ///< \f$A_{1}\f$ PSS numerator coefficient (first block)
        A2,     ///< \f$A_{2}\f$ PSS numerator coefficient (first block)
        A3,     ///< \f$A_{3}\f$ PSS denominator coefficient (first block)
        A4,     ///< \f$A_{4}\f$ PSS denominator coefficient (first block)
        A5,     ///< \f$A_{5}\f$ PSS numerator coefficient (second block)
        A6,     ///< \f$A_{6}\f$ PSS numerator coefficient (second block)
        A7,     ///< \f$A_{7}\f$ PSS denominator coefficient (second block)
        A8,     ///< \f$A_{8}\f$ PSS denominator coefficient (second block)
        VSTMAX, ///< \f$V_{\mathrm{STMAX}}\f$ Maximum PSS output
        VSTMIN, ///< \f$V_{\mathrm{STMIN}}\f$ Minimum PSS output
        PPSSon, ///< \f$P_{\mathrm{PSSon}}\f$ Generator MW threshold for PSS activation
        PPSSoff ///< \f$P_{\mathrm{PSSoff}}\f$ Generator MW threshold for PSS de-activation
      };

      /// Bus keys for the PSS3C model; deferred until port integration.
      enum class Pss3cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the PSS3C model; deferred until port integration.
      enum class Pss3cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the PSS3C model; deferred until port integration.
      enum class Pss3cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the PSS3C model; deferred until implementation.
      enum class Pss3cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for PSS3C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Pss3cData : public ComponentData<real_type,
                                              index_type,
                                              Pss3cParameters,
                                              Pss3cBuses,
                                              Pss3cSignalInputs,
                                              Pss3cSignalOutputs,
                                              Pss3cMonitorableVariables>
      {
        Pss3cData() = default;

        using Parameters           = Pss3cParameters;
        using Buses                = Pss3cBuses;
        using SignalInputs         = Pss3cSignalInputs;
        using SignalOutputs        = Pss3cSignalOutputs;
        using MonitorableVariables = Pss3cMonitorableVariables;
      };
    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit
