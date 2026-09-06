/**
 * @file IeeestData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for IEEEST Controller
 */

#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      /**
       * @brief Parameter keys for IEEEST Stabilizer model.
       */
      enum class IeeestParameters
      {
        A1,     ///< \f$A_1\f$ Notch-filter denominator coefficient
        A2,     ///< \f$A_2\f$ Notch-filter denominator coefficient
        A3,     ///< \f$A_3\f$ Notch-filter denominator coefficient
        A4,     ///< \f$A_4\f$ Notch-filter denominator coefficient
        A5,     ///< \f$A_5\f$ Notch-filter numerator coefficient
        A6,     ///< \f$A_6\f$ Notch-filter numerator coefficient
        T1,     ///< \f$T_1\f$ Lead-lag 1 numerator time constant
        T2,     ///< \f$T_2\f$ Lead-lag 1 denominator time constant
        T3,     ///< \f$T_3\f$ Lead-lag 2 numerator time constant
        T4,     ///< \f$T_4\f$ Lead-lag 2 denominator time constant
        T5,     ///< \f$T_5\f$ Washout numerator time constant
        T6,     ///< \f$T_6\f$ Washout denominator time constant
        Ks,     ///< \f$K_s\f$ Stabilizer gain
        Lsmin,  ///< \f$L_s^{\min}\f$ Minimum stabilizer output limit
        Lsmax,  ///< \f$L_s^{\max}\f$ Maximum stabilizer output limit
        Vcl,    ///< \f$V_{\mathrm{cl}}\f$ Lower voltage-cutout threshold; zero disables
        Vcu,    ///< \f$V_{\mathrm{cu}}\f$ Upper voltage-cutout threshold; zero disables
        Tdelay, ///< \f$T_{\mathrm{delay}}\f$ PSLF delay extension; only zero is supported
      };

      /**
       * @brief Signal input keys for IEEEST Stabilizer model.
       */
      enum class IeeestInputs : size_t
      {
        input, ///< Generic stabilizer input signal, alternative to speed
        speed, ///< Absolute per-unit rotor speed; model subtracts one
        vct,   ///< Compensated voltage magnitude in per unit, required when cutout is enabled
        SIZE
      };

      /**
       * @brief Signal output keys for IEEEST Stabilizer model.
       */
      enum class IeeestOutputs : size_t
      {
        output, ///< Unique ID of the stabilizer output signal
        SIZE
      };

      /**
       * @brief Monitorable variables for IEEEST Stabilizer model.
       */
      enum class IeeestMonitorableVariables
      {
        vss, ///< \f$V_{\mathrm{ss}}\f$ Limited stabilizer signal and model output
      };

      /**
       * @brief Contains modeling data for a IEEEST Stabilizer model.
       *
       * @tparam real_type  Real parameter data type
       * @tparam index_type Integer parameter data type
       */
      template <typename real_type, typename index_type>
      struct IeeestData : public ComponentData<real_type,
                                               index_type,
                                               IeeestParameters,
                                               IeeestInputs,
                                               IeeestOutputs,
                                               IeeestMonitorableVariables>
      {
        IeeestData() = default;

        using Parameters           = IeeestParameters;
        using Inputs               = IeeestInputs;
        using Outputs              = IeeestOutputs;
        using MonitorableVariables = IeeestMonitorableVariables;
      };

    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
