/**
 * @file Uel2cData.hpp
 * @brief Modeling data for the UEL2C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Limiter
    {
      /// Parameter keys for the UEL2C model.
      enum class Uel2cParameters
      {
        TUP,      ///< \f$T_{\mathrm{UP}}\f$ UEL real power filter time constant
        TUQ,      ///< \f$T_{\mathrm{UQ}}\f$ UEL reactive power filter time constant
        TUV,      ///< \f$T_{\mathrm{UV}}\f$ UEL voltage filter time constant
        Vbias,    ///< \f$V_{\mathrm{bias}}\f$ UEL voltage bias
        K1,       ///< \f$K_{1}\f$ Voltage exponent for real power input to UEL table
        K2,       ///< \f$K_{2}\f$ Voltage exponent for reactive power output of UEL table
        KUF,      ///< \f$K_{\mathrm{UF}}\f$ UEL excitation system stabilizer gain
        TQref,    ///< \f$T_{\mathrm{Qref}}\f$ UEL reactive power reference time constant
        Kfix,     ///< \f$K_{\mathrm{fix}}\f$ UEL fixed gain reduction factor
        Tadj,     ///< \f$T_{\mathrm{adj}}\f$ UEL adjustable gain reduction time constant
        SW1,      ///< \f$\mathrm{SW}_{1}\f$ UEL logic switch for adjustable gain reduction
        KUI,      ///< \f$K_{\mathrm{UI}}\f$ UEL integral gain
        KUL,      ///< \f$K_{\mathrm{UL}}\f$ UEL proportional gain
        VUImax,   ///< \f$V_{\mathrm{UImax}}\f$ UEL PI control maximum output
        VUImin,   ///< \f$V_{\mathrm{UImin}}\f$ UEL PI control minimum output
        TU1,      ///< \f$T_{\mathrm{U1}}\f$ UEL numerator (lead) time constant (first block)
        TU2,      ///< \f$T_{\mathrm{U2}}\f$ UEL denominator (lag) time constant (first block)
        TU3,      ///< \f$T_{\mathrm{U3}}\f$ UEL numerator (lead) time constant (second block)
        TU4,      ///< \f$T_{\mathrm{U4}}\f$ UEL denominator (lag) time constant (second block)
        VUELmax1, ///< \f$V_{\mathrm{UELmax1}}\f$ UEL maximum output
        VUELmin1, ///< \f$V_{\mathrm{UELmin1}}\f$ UEL minimum output
        VUELmax2, ///< \f$V_{\mathrm{UELmax2}}\f$ UEL maximum output
        VUELmin2, ///< \f$V_{\mathrm{UELmin2}}\f$ UEL minimum output
        P0,       ///< \f$P_{0}\f$ UEL lookup table real power (first point)
        Q0,       ///< \f$Q_{0}\f$ UEL lookup table reactive power (first point)
        P1,       ///< \f$P_{1}\f$ UEL lookup table real power (second point)
        Q1,       ///< \f$Q_{1}\f$ UEL lookup table reactive power (second point)
        P2,       ///< \f$P_{2}\f$ UEL lookup table real power (third point)
        Q2,       ///< \f$Q_{2}\f$ UEL lookup table reactive power (third point)
        P3,       ///< \f$P_{3}\f$ UEL lookup table real power (fourth point)
        Q3,       ///< \f$Q_{3}\f$ UEL lookup table reactive power (fourth point)
        P4,       ///< \f$P_{4}\f$ UEL lookup table real power (fifth point)
        Q4        ///< \f$Q_{4}\f$ UEL lookup table reactive power (fifth point)
      };

      /// Bus keys for the UEL2C model; deferred until port integration.
      enum class Uel2cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the UEL2C model; deferred until port integration.
      enum class Uel2cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the UEL2C model; deferred until port integration.
      enum class Uel2cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the UEL2C model; deferred until implementation.
      enum class Uel2cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for UEL2C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Uel2cData : public ComponentData<real_type,
                                              index_type,
                                              Uel2cParameters,
                                              Uel2cBuses,
                                              Uel2cSignalInputs,
                                              Uel2cSignalOutputs,
                                              Uel2cMonitorableVariables>
      {
        Uel2cData() = default;

        using Parameters           = Uel2cParameters;
        using Buses                = Uel2cBuses;
        using SignalInputs         = Uel2cSignalInputs;
        using SignalOutputs        = Uel2cSignalOutputs;
        using MonitorableVariables = Uel2cMonitorableVariables;
      };
    } // namespace Limiter
  } // namespace PhasorDynamics
} // namespace GridKit
