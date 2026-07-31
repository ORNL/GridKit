/**
 * @file Pss7cData.hpp
 * @brief Modeling data for the PSS7C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Stabilizer
    {
      /// Parameter keys for the PSS7C model.
      enum class Pss7cParameters
      {
        KS1,     ///< \f$K_{\mathrm{S1}}\f$ Power system stabilizer main gain
        KS2,     ///< \f$K_{\mathrm{S2}}\f$ Power system stabilizer gain
        KS3,     ///< \f$K_{\mathrm{S3}}\f$ Power system stabilizer gain
        T6,      ///< \f$T_{6}\f$ PSS transducer time constant
        T7,      ///< \f$T_{7}\f$ PSS transducer time constant
        Tw1,     ///< \f$T_{\mathrm{w1}}\f$ PSS washout time constant
        Tw2,     ///< \f$T_{\mathrm{w2}}\f$ PSS washout time constant
        Tw3,     ///< \f$T_{\mathrm{w3}}\f$ PSS washout time constant
        Tw4,     ///< \f$T_{\mathrm{w4}}\f$ PSS washout time constant
        T8,      ///< \f$T_{8}\f$ PSS transducer time constant
        T9,      ///< \f$T_{9}\f$ PSS washout time constant
        M,       ///< \f$M\f$ Denominator exponent for ramp-track filter
        N,       ///< \f$N\f$ Overall exponent for ramp-track filter
        K0,      ///< \f$K_{0}\f$ PSS canonical gain 0
        K1,      ///< \f$K_{1}\f$ PSS canonical gain 1
        K2,      ///< \f$K_{2}\f$ PSS canonical gain 2
        K3,      ///< \f$K_{3}\f$ PSS canonical gain 3
        K4,      ///< \f$K_{4}\f$ PSS canonical gain 4
        Ki3,     ///< \f$K_{\mathrm{i3}}\f$ PSS third block gain
        Ki4,     ///< \f$K_{\mathrm{i4}}\f$ PSS fourth block gain
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

      /// Bus keys for the PSS7C model; deferred until port integration.
      enum class Pss7cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the PSS7C model; deferred until port integration.
      enum class Pss7cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the PSS7C model; deferred until port integration.
      enum class Pss7cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the PSS7C model; deferred until implementation.
      enum class Pss7cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for PSS7C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Pss7cData : public ComponentData<real_type,
                                              index_type,
                                              Pss7cParameters,
                                              Pss7cBuses,
                                              Pss7cSignalInputs,
                                              Pss7cSignalOutputs,
                                              Pss7cMonitorableVariables>
      {
        Pss7cData() = default;

        using Parameters           = Pss7cParameters;
        using Buses                = Pss7cBuses;
        using SignalInputs         = Pss7cSignalInputs;
        using SignalOutputs        = Pss7cSignalOutputs;
        using MonitorableVariables = Pss7cMonitorableVariables;
      };
    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit
