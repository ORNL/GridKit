/**
 * @file Pss6cData.hpp
 * @brief Modeling data for the PSS6C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Stabilizer
    {
      /// Parameter keys for the PSS6C model.
      enum class Pss6cParameters
      {
        KS1,     ///< \f$K_{\mathrm{S1}}\f$ PSS gain for input channel 1
        T1,      ///< \f$T_{1}\f$ PSS transducer time constant for input channel 1
        T3,      ///< \f$T_{3}\f$ PSS time constant for input channel 1
        KS2,     ///< \f$K_{\mathrm{S2}}\f$ PSS gain for input channel 2
        Macc,    ///< \f$M_{\mathrm{acc}}\f$ PSS washout time constant for input channel 2
        T2,      ///< \f$T_{2}\f$ PSS transducer time constant for input channel 2
        T4,      ///< \f$T_{4}\f$ PSS time constant for input channel 2
        TD,      ///< \f$T_{\mathrm{D}}\f$ PSS washout time constant
        K0,      ///< \f$K_{0}\f$ PSS canonical gain 0
        K1,      ///< \f$K_{1}\f$ PSS canonical gain 1
        K2,      ///< \f$K_{2}\f$ PSS canonical gain 2
        K3,      ///< \f$K_{3}\f$ PSS canonical gain 3
        K4,      ///< \f$K_{4}\f$ PSS canonical gain 4
        Ki3,     ///< \f$K_{\mathrm{i3}}\f$ PSS third block gain
        Ki4,     ///< \f$K_{\mathrm{i4}}\f$ PSS fourth block gain
        Ks,      ///< \f$K_{\mathrm{s}}\f$ PSS main gain
        Ti1,     ///< \f$T_{\mathrm{i1}}\f$ PSS time constant (first block)
        Ti2,     ///< \f$T_{\mathrm{i2}}\f$ PSS time constant (second block)
        Ti3,     ///< \f$T_{\mathrm{i3}}\f$ PSS time constant (third block)
        Ti4,     ///< \f$T_{\mathrm{i4}}\f$ PSS time constant (fourth block)
        VSI1max, ///< \f$V_{\mathrm{SI1max}}\f$ Input signal #1 maximum limit
        VSI1min, ///< \f$V_{\mathrm{SI1min}}\f$ Input signal #1 minimum limit
        VSI2max, ///< \f$V_{\mathrm{SI2max}}\f$ Input signal #2 maximum limit
        VSI2min, ///< \f$V_{\mathrm{SI2min}}\f$ Input signal #2 minimum limit
        VSTMAX,  ///< \f$V_{\mathrm{STMAX}}\f$ Maximum PSS output
        VSTMIN,  ///< \f$V_{\mathrm{STMIN}}\f$ Minimum PSS output
        PPSSon,  ///< \f$P_{\mathrm{PSSon}}\f$ Generator MW threshold for PSS activation
        PPSSoff  ///< \f$P_{\mathrm{PSSoff}}\f$ Generator MW threshold for PSS de-activation
      };

      /// Bus keys for the PSS6C model; deferred until port integration.
      enum class Pss6cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the PSS6C model; deferred until port integration.
      enum class Pss6cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the PSS6C model; deferred until port integration.
      enum class Pss6cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the PSS6C model; deferred until implementation.
      enum class Pss6cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for PSS6C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Pss6cData : public ComponentData<real_type,
                                              index_type,
                                              Pss6cParameters,
                                              Pss6cBuses,
                                              Pss6cSignalInputs,
                                              Pss6cSignalOutputs,
                                              Pss6cMonitorableVariables>
      {
        Pss6cData() = default;

        using Parameters           = Pss6cParameters;
        using Buses                = Pss6cBuses;
        using SignalInputs         = Pss6cSignalInputs;
        using SignalOutputs        = Pss6cSignalOutputs;
        using MonitorableVariables = Pss6cMonitorableVariables;
      };
    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit
