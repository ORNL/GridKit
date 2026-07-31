/**
 * @file Pss4cData.hpp
 * @brief Modeling data for the PSS4C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Stabilizer
    {
      /// Parameter keys for the PSS4C model.
      enum class Pss4cParameters
      {
        KVL,    ///< \f$K_{\mathrm{VL}}\f$ Very low band gain
        KVL1,   ///< \f$K_{\mathrm{VL1}}\f$ Very low band differential filter gain
        KVL11,  ///< \f$K_{\mathrm{VL11}}\f$ Very low band first lead-lag block coefficient
        TVL1,   ///< \f$T_{\mathrm{VL1}}\f$ Very low band numerator time constant (first lead-lag block)
        TVL2,   ///< \f$T_{\mathrm{VL2}}\f$ Very low band denominator time constant (first lead-lag block)
        TVL3,   ///< \f$T_{\mathrm{VL3}}\f$ Very low band numerator time constant (second lead-lag block)
        TVL4,   ///< \f$T_{\mathrm{VL4}}\f$ Very low band denominator time constant (second lead-lag block)
        TVL5,   ///< \f$T_{\mathrm{VL5}}\f$ Very low band numerator time constant (third lead-lag block)
        TVL6,   ///< \f$T_{\mathrm{VL6}}\f$ Very low band denominator time constant (third lead-lag block)
        KVL2,   ///< \f$K_{\mathrm{VL2}}\f$ Very low band differential filter gain
        KVL17,  ///< \f$K_{\mathrm{VL17}}\f$ Very low band first lead-lag block coefficient
        TVL7,   ///< \f$T_{\mathrm{VL7}}\f$ Very low band numerator time constant (first lead-lag block)
        TVL8,   ///< \f$T_{\mathrm{VL8}}\f$ Very low band denominator time constant (first lead-lag block)
        TVL9,   ///< \f$T_{\mathrm{VL9}}\f$ Very low band numerator time constant (second lead-lag block)
        TVL10,  ///< \f$T_{\mathrm{VL10}}\f$ Very low band denominator time constant (second lead-lag block)
        TVL11,  ///< \f$T_{\mathrm{VL11}}\f$ Very low band numerator time constant (third lead-lag block)
        TVL12,  ///< \f$T_{\mathrm{VL12}}\f$ Very low band denominator time constant (third lead-lag block)
        VVLmax, ///< \f$V_{\mathrm{VLmax}}\f$ Very low band upper limit
        VVLmin, ///< \f$V_{\mathrm{VLmin}}\f$ Very low band lower limit
        KL,     ///< \f$K_{\mathrm{L}}\f$ Low band gain
        KL1,    ///< \f$K_{\mathrm{L1}}\f$ Low band differential filter gain
        KL11,   ///< \f$K_{\mathrm{L11}}\f$ Low band first lead-lag block coefficient
        TL1,    ///< \f$T_{\mathrm{L1}}\f$ Low band numerator time constant (first lead-lag block)
        TL2,    ///< \f$T_{\mathrm{L2}}\f$ Low band denominator time constant (first lead-lag block)
        TL3,    ///< \f$T_{\mathrm{L3}}\f$ Low band numerator time constant (second lead-lag block)
        TL4,    ///< \f$T_{\mathrm{L4}}\f$ Low band denominator time constant (second lead-lag block)
        TL5,    ///< \f$T_{\mathrm{L5}}\f$ Low band numerator time constant (third lead-lag block)
        TL6,    ///< \f$T_{\mathrm{L6}}\f$ Low band denominator time constant (third lead-lag block)
        KL2,    ///< \f$K_{\mathrm{L2}}\f$ Low band differential filter gain
        KL17,   ///< \f$K_{\mathrm{L17}}\f$ Low band first lead-lag block coefficient
        TL7,    ///< \f$T_{\mathrm{L7}}\f$ Low band numerator time constant (first lead-lag block)
        TL8,    ///< \f$T_{\mathrm{L8}}\f$ Low band denominator time constant (first lead-lag block)
        TL9,    ///< \f$T_{\mathrm{L9}}\f$ Low band numerator time constant (second lead-lag block)
        TL10,   ///< \f$T_{\mathrm{L10}}\f$ Low band denominator time constant (second lead-lag block)
        TL11,   ///< \f$T_{\mathrm{L11}}\f$ Low band numerator time constant (third lead-lag block)
        TL12,   ///< \f$T_{\mathrm{L12}}\f$ Low band denominator time constant (third lead-lag block)
        VLmax,  ///< \f$V_{\mathrm{Lmax}}\f$ Low band upper limit
        VLmin,  ///< \f$V_{\mathrm{Lmin}}\f$ Low band lower limit
        KI,     ///< \f$K_{\mathrm{I}}\f$ Intermediate band gain
        KI1,    ///< \f$K_{\mathrm{I1}}\f$ Intermediate band differential filter gain
        KI11,   ///< \f$K_{\mathrm{I11}}\f$ Intermediate band first lead-lag block coefficient
        TI1,    ///< \f$T_{\mathrm{I1}}\f$ Intermediate band numerator time constant (first lead-lag block)
        TI2,    ///< \f$T_{\mathrm{I2}}\f$ Intermediate band denominator time constant (first block)
        TI3,    ///< \f$T_{\mathrm{I3}}\f$ Intermediate band numerator time constant (second lead-lag block)
        TI4,    ///< \f$T_{\mathrm{I4}}\f$ Intermediate band denominator time constant (second block)
        TI5,    ///< \f$T_{\mathrm{I5}}\f$ Intermediate band numerator time constant (third lead-lag block)
        TI6,    ///< \f$T_{\mathrm{I6}}\f$ Intermediate band denominator time constant (third block)
        KI2,    ///< \f$K_{\mathrm{I2}}\f$ Intermediate band differential filter gain
        KI17,   ///< \f$K_{\mathrm{I17}}\f$ Intermediate band first lead-lag block coefficient
        TI7,    ///< \f$T_{\mathrm{I7}}\f$ Intermediate band numerator time constant (first lead-lag block)
        TI8,    ///< \f$T_{\mathrm{I8}}\f$ Intermediate band denominator time constant (first block)
        TI9,    ///< \f$T_{\mathrm{I9}}\f$ Intermediate band numerator time constant (second lead-lag block)
        TI10,   ///< \f$T_{\mathrm{I10}}\f$ Intermediate band denominator time constant (second block)
        TI11,   ///< \f$T_{\mathrm{I11}}\f$ Intermediate band numerator time constant (third lead-lag block)
        TI12,   ///< \f$T_{\mathrm{I12}}\f$ Intermediate band denominator time constant (third block)
        VImax,  ///< \f$V_{\mathrm{Imax}}\f$ Intermediate band upper limit
        VImin,  ///< \f$V_{\mathrm{Imin}}\f$ Intermediate band lower limit
        KH,     ///< \f$K_{\mathrm{H}}\f$ High band gain
        KH1,    ///< \f$K_{\mathrm{H1}}\f$ High band differential filter gain
        KH11,   ///< \f$K_{\mathrm{H11}}\f$ High band first lead-lag block coefficient
        TH1,    ///< \f$T_{\mathrm{H1}}\f$ High band numerator time constant (first lead-lag block)
        TH2,    ///< \f$T_{\mathrm{H2}}\f$ High band denominator time constant (first lead-lag block)
        TH3,    ///< \f$T_{\mathrm{H3}}\f$ High band numerator time constant (second lead-lag block)
        TH4,    ///< \f$T_{\mathrm{H4}}\f$ High band denominator time constant (second lead-lag block)
        TH5,    ///< \f$T_{\mathrm{H5}}\f$ High band numerator time constant (third lead-lag block)
        TH6,    ///< \f$T_{\mathrm{H6}}\f$ High band denominator time constant (third lead-lag block)
        KH2,    ///< \f$K_{\mathrm{H2}}\f$ High band differential filter gain
        KH17,   ///< \f$K_{\mathrm{H17}}\f$ High band first lead-lag block coefficient
        TH7,    ///< \f$T_{\mathrm{H7}}\f$ High band numerator time constant (first lead-lag block)
        TH8,    ///< \f$T_{\mathrm{H8}}\f$ High band denominator time constant (first lead-lag block)
        TH9,    ///< \f$T_{\mathrm{H9}}\f$ High band numerator time constant (second lead-lag block)
        TH10,   ///< \f$T_{\mathrm{H10}}\f$ High band denominator time constant (second lead-lag block)
        TH11,   ///< \f$T_{\mathrm{H11}}\f$ High band numerator time constant (third lead-lag block)
        TH12,   ///< \f$T_{\mathrm{H12}}\f$ High band denominator time constant (third lead-lag block)
        VHmax,  ///< \f$V_{\mathrm{Hmax}}\f$ High band upper limit
        VHmin,  ///< \f$V_{\mathrm{Hmin}}\f$ High band lower limit
        VSTMAX, ///< \f$V_{\mathrm{STMAX}}\f$ Maximum PSS output
        VSTMIN  ///< \f$V_{\mathrm{STMIN}}\f$ Minimum PSS output
      };

      /// Bus keys for the PSS4C model; deferred until port integration.
      enum class Pss4cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the PSS4C model; deferred until port integration.
      enum class Pss4cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the PSS4C model; deferred until port integration.
      enum class Pss4cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the PSS4C model; deferred until implementation.
      enum class Pss4cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for PSS4C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Pss4cData : public ComponentData<real_type,
                                              index_type,
                                              Pss4cParameters,
                                              Pss4cBuses,
                                              Pss4cSignalInputs,
                                              Pss4cSignalOutputs,
                                              Pss4cMonitorableVariables>
      {
        Pss4cData() = default;

        using Parameters           = Pss4cParameters;
        using Buses                = Pss4cBuses;
        using SignalInputs         = Pss4cSignalInputs;
        using SignalOutputs        = Pss4cSignalOutputs;
        using MonitorableVariables = Pss4cMonitorableVariables;
      };
    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit
