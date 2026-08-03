/**
 * @file IeeestData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for IEEEST Stabilizer
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Stabilizer
    {
      /**
       * @brief Parameter keys for IEEEST Stabilizer model.
       */
      enum class IeeestParameters
      {
        A1,     ///< \f$A_1\f$ Notch filter denominator coefficient [sec]
        A2,     ///< \f$A_2\f$ Notch filter denominator coefficient [sec^2]
        A3,     ///< \f$A_3\f$ Notch filter denominator coefficient [sec]
        A4,     ///< \f$A_4\f$ Notch filter denominator coefficient [sec^2]
        A5,     ///< \f$A_5\f$ Notch filter numerator coefficient [sec]
        A6,     ///< \f$A_6\f$ Notch filter numerator coefficient [sec^2]
        T1,     ///< \f$T_1\f$ Lead-lag 1 numerator time constant [sec]
        T2,     ///< \f$T_2\f$ Lead-lag 1 denominator time constant [sec]
        T3,     ///< \f$T_3\f$ Lead-lag 2 numerator time constant [sec]
        T4,     ///< \f$T_4\f$ Lead-lag 2 denominator time constant [sec]
        T5,     ///< \f$T_5\f$ Washout numerator time constant [sec]
        T6,     ///< \f$T_6\f$ Washout denominator time constant [sec]
        Ks,     ///< \f$K_s\f$ Stabilizer gain [p.u.]
        Lsmin,  ///< \f$L_s^{\min}\f$ Minimum stabilizer output limit [p.u.]
        Lsmax,  ///< \f$L_s^{\max}\f$ Maximum stabilizer output limit [p.u.]
        Vcl,    ///< \f$V_{cl}\f$ Lower input cutout threshold [p.u.], accepted and unused
        Vcu,    ///< \f$V_{cu}\f$ Upper input cutout threshold [p.u.], accepted and unused
        Tdelay, ///< \f$T_{delay}\f$ Input time delay [sec], accepted and unused
      };

      /// Buses for the IEEEST Stabilizer model. The stabilizer has no
      /// electrical terminal.
      enum class IeeestBuses : size_t
      {
        SIZE
      };

      /// Signal inputs for the IEEEST Stabilizer model.
      enum class IeeestSignalInputs : size_t
      {
        input, ///< \f$u\f$ Required stabilizer input signal ID
        SIZE
      };

      /// Signal outputs for the IEEEST Stabilizer model.
      enum class IeeestSignalOutputs : size_t
      {
        output, ///< \f$V_{ss}\f$ Optional stabilizer output signal ID
        SIZE
      };

      /// Variables available through the monitor interface.
      enum class IeeestMonitorableVariables
      {
        vss, ///< \f$V_{ss}\f$ Stabilizer output after the limiter [p.u.]
      };

      /**
       * @brief Model data for the IEEEST Stabilizer: parameter values, the
       *        input and output signals, and monitored-variable selections.
       *
       * @tparam real_type  Real parameter data type
       * @tparam index_type Integer parameter data type
       *
       * @see Ieeest
       */
      template <typename real_type, typename index_type>
      struct IeeestData : public ComponentData<real_type,
                                               index_type,
                                               IeeestParameters,
                                               IeeestBuses,
                                               IeeestSignalInputs,
                                               IeeestSignalOutputs,
                                               IeeestMonitorableVariables>
      {
        IeeestData() = default;

        using Parameters           = IeeestParameters;
        using Buses                = IeeestBuses;
        using SignalInputs         = IeeestSignalInputs;
        using SignalOutputs        = IeeestSignalOutputs;
        using MonitorableVariables = IeeestMonitorableVariables;
      };

    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit
