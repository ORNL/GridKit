/**
 * @file Pss2cData.hpp
 * @brief Modeling data for the PSS2C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Stabilizer
    {
      /// Parameter keys for the PSS2C model.
      enum class Pss2cParameters
      {
        KS1,     ///< \f$K_{\mathrm{S1}}\f$ Power system stabilizer gain
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
        M,       ///< \f$M\f$ PSS transducer time constant
        N,       ///< \f$N\f$ PSS washout time constant
        T1,      ///< \f$T_{1}\f$ PSS numerator (lead) compensating time constant (first block)
        T2,      ///< \f$T_{2}\f$ PSS denominator (lag) compensating time constant (first block)
        T3,      ///< \f$T_{3}\f$ PSS numerator (lead) compensating time constant (second block)
        T4,      ///< \f$T_{4}\f$ PSS denominator (lag) compensating time constant (second block)
        T10,     ///< \f$T_{10}\f$ PSS numerator (lead) compensating time constant (third block)
        T11,     ///< \f$T_{11}\f$ PSS denominator (lag) compensating time constant (third block)
        T12,     ///< \f$T_{12}\f$ PSS numerator (lead) compensating time constant (fourth block)
        T13,     ///< \f$T_{13}\f$ PSS denominator (lag) compensating time constant (fourth block)
        VSTmax,  ///< \f$V_{\mathrm{STmax}}\f$ Maximum PSS output
        VSTmin,  ///< \f$V_{\mathrm{STmin}}\f$ Minimum PSS output
        VSI1max, ///< \f$V_{\mathrm{SI1max}}\f$ Input signal #1 maximum limit
        VSI1min, ///< \f$V_{\mathrm{SI1min}}\f$ Input signal #1 minimum limit
        VSI2max, ///< \f$V_{\mathrm{SI2max}}\f$ Input signal #2 maximum limit
        VSI2min, ///< \f$V_{\mathrm{SI2min}}\f$ Input signal #2 minimum limit
        PPSSon,  ///< \f$P_{\mathrm{PSSon}}\f$ Generator MW threshold for PSS activation
        PPSSoff  ///< \f$P_{\mathrm{PSSoff}}\f$ Generator MW threshold for PSS de-activation
      };

      /// Bus keys for the PSS2C model; deferred until port integration.
      enum class Pss2cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the PSS2C model; deferred until port integration.
      enum class Pss2cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the PSS2C model; deferred until port integration.
      enum class Pss2cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the PSS2C model; deferred until implementation.
      enum class Pss2cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for PSS2C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Pss2cData : public ComponentData<real_type,
                                              index_type,
                                              Pss2cParameters,
                                              Pss2cBuses,
                                              Pss2cSignalInputs,
                                              Pss2cSignalOutputs,
                                              Pss2cMonitorableVariables>
      {
        Pss2cData() = default;

        using Parameters           = Pss2cParameters;
        using Buses                = Pss2cBuses;
        using SignalInputs         = Pss2cSignalInputs;
        using SignalOutputs        = Pss2cSignalOutputs;
        using MonitorableVariables = Pss2cMonitorableVariables;
      };
    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit
