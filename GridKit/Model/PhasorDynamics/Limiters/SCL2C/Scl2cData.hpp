/**
 * @file Scl2cData.hpp
 * @brief Modeling data for the SCL2C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Limiter
    {
      /// Parameter keys for the SCL2C model.
      enum class Scl2cParameters
      {
        TB1oel,   ///< \f$T_{\mathrm{B1oel}}\f$ Overexcited regulator denominator (lag) time constant 1
        TC1oel,   ///< \f$T_{\mathrm{C1oel}}\f$ Overexcited regulator numerator (lead) time constant 1
        TB2oel,   ///< \f$T_{\mathrm{B2oel}}\f$ Overexcited regulator denominator (lag) time constant 2
        TC2oel,   ///< \f$T_{\mathrm{C2oel}}\f$ Overexcited regulator numerator (lead) time constant 2
        KPoel,    ///< \f$K_{\mathrm{Poel}}\f$ Overexcited PID regulator proportional gain
        KIoel,    ///< \f$K_{\mathrm{Ioel}}\f$ Overexcited PID regulator integral gain
        KDoel,    ///< \f$K_{\mathrm{Doel}}\f$ Overexcited PID regulator differential gain
        TDoel,    ///< \f$T_{\mathrm{Doel}}\f$ Overexcited PID regulator differential time constant
        VOELmax3, ///< \f$V_{\mathrm{OELmax3}}\f$ Maximum OEL PID output limit
        VOELmin3, ///< \f$V_{\mathrm{OELmin3}}\f$ Minimum OEL PID output limit
        VOELmax2, ///< \f$V_{\mathrm{OELmax2}}\f$ Maximum OEL lead-lag 1 output limit
        VOELmin2, ///< \f$V_{\mathrm{OELmin2}}\f$ Minimum OEL lead-lag 1 output limit
        VOELmax1, ///< \f$V_{\mathrm{OELmax1}}\f$ Maximum OEL output limit
        VOELmin1, ///< \f$V_{\mathrm{OELmin1}}\f$ Minimum OEL output limit
        TB1uel,   ///< \f$T_{\mathrm{B1uel}}\f$ Underexcited regulator denominator (lag) time constant 1
        TC1uel,   ///< \f$T_{\mathrm{C1uel}}\f$ Underexcited regulator numerator (lead) time constant 1
        TB2uel,   ///< \f$T_{\mathrm{B2uel}}\f$ Underexcited regulator denominator (lag) time constant 2
        TC2uel,   ///< \f$T_{\mathrm{C2uel}}\f$ Underexcited regulator numerator (lead) time constant 2
        KPuel,    ///< \f$K_{\mathrm{Puel}}\f$ Underexcited PID regulator proportional gain
        KIuel,    ///< \f$K_{\mathrm{Iuel}}\f$ Underexcited PID regulator integral gain
        KDuel,    ///< \f$K_{\mathrm{Duel}}\f$ Underexcited PID regulator differential gain
        TDuel,    ///< \f$T_{\mathrm{Duel}}\f$ Underexcited PID regulator differential time constant
        VUELmax3, ///< \f$V_{\mathrm{UELmax3}}\f$ Maximum UEL PID output limit
        VUELmin3, ///< \f$V_{\mathrm{UELmin3}}\f$ Minimum UEL PID output limit
        VUELmax2, ///< \f$V_{\mathrm{UELmax2}}\f$ Maximum UEL lead-lag 1 output limit
        VUELmin2, ///< \f$V_{\mathrm{UELmin2}}\f$ Minimum UEL lead-lag 1 output limit
        VUELmax1, ///< \f$V_{\mathrm{UELmax1}}\f$ Maximum UEL output limit
        VUELmin1, ///< \f$V_{\mathrm{UELmin1}}\f$ Minimum UEL output limit
        Ireset,   ///< \f$I_{\mathrm{reset}}\f$ SCL reset-reference, if inactive
        TenOEL,   ///< \f$T_{\mathrm{enOEL}}\f$ Overexcited activation delay time
        TenUEL,   ///< \f$T_{\mathrm{enUEL}}\f$ Underexcited activation delay time
        Toff,     ///< \f$T_{\mathrm{off}}\f$ SCL reset delay time
        ITHoff,   ///< \f$I_{\mathrm{THoff}}\f$ SCL reset threshold value
        TIQoel,   ///< \f$T_{\mathrm{IQoel}}\f$ Overexcited reactive current time constant
        KIQoel,   ///< \f$K_{\mathrm{IQoel}}\f$ Overexcited reactive current scaling factor
        TIPoel,   ///< \f$T_{\mathrm{IPoel}}\f$ Overexcited active current time constant
        KIPoel,   ///< \f$K_{\mathrm{IPoel}}\f$ Overexcited active current scaling factor
        TIQuel,   ///< \f$T_{\mathrm{IQuel}}\f$ Underexcited reactive current time constant
        KIQuel,   ///< \f$K_{\mathrm{IQuel}}\f$ Underexcited reactive current scaling factor
        TIPuel,   ///< \f$T_{\mathrm{IPuel}}\f$ Underexcited active current time constant
        KIPuel,   ///< \f$K_{\mathrm{IPuel}}\f$ Underexcited active current scaling factor
        TITscl,   ///< \f$T_{\mathrm{ITscl}}\f$ Stator current transducer time constant
        ITFpu,    ///< \f$I_{\mathrm{TFpu}}\f$ SCL thermal reference for inverse time calculations
        Iinst,    ///< \f$I_{\mathrm{inst}}\f$ SCL instantaneous stator current limit
        IinstUEL, ///< \f$I_{\mathrm{instUEL}}\f$ Underexcited region instantaneous stator current limit
        Ilim,     ///< \f$I_{\mathrm{lim}}\f$ SCL thermal stator current limit
        TAoel,    ///< \f$T_{\mathrm{Aoel}}\f$ SCL reference filter time constant
        c1,       ///< \f$c_{1}\f$ SCL exponent for calculation of I ERRinv1
        K1,       ///< \f$K_{1}\f$ SCL gain for calculation of I ERRinv1
        c2,       ///< \f$c_{2}\f$ SCL exponent for calculation of I c ERRinv2
        K2,       ///< \f$K_{2}\f$ SCL gain for calculation of I c ERRinv2
        VINVmax,  ///< \f$V_{\mathrm{INVmax}}\f$ SCL maximum inverse time output
        VINVmin,  ///< \f$V_{\mathrm{INVmin}}\f$ SCL minimum inverse time output
        Fixedru,  ///< \f$\mathrm{Fixed}_{\mathrm{ru}}\f$ SCL fixed delay time output
        Fixedrd,  ///< \f$\mathrm{Fixed}_{\mathrm{rd}}\f$ SCL fixed cooling-down time output
        TSCL,     ///< \f$T_{\mathrm{SCL}}\f$ SCL timer reference
        Tmax,     ///< \f$T_{\mathrm{max}}\f$ SCL timer maximum level
        Tmin,     ///< \f$T_{\mathrm{min}}\f$ SCL timer minimum level
        KFB,      ///< \f$K_{\mathrm{FB}}\f$ SCL timer feedback gain
        SW1,      ///< \f$\mathrm{SW1}\f$ OEL reference ramp logic selection
        Krd,      ///< \f$K_{\mathrm{rd}}\f$ SCL reference ramp-down rate
        Kru,      ///< \f$K_{\mathrm{ru}}\f$ SCL reference ramp-up rate
        KZRU,     ///< \f$K_{\mathrm{ZRU}}\f$ SCL thermal reference release threshold
        TVTscl,   ///< \f$T_{\mathrm{VTscl}}\f$ Terminal voltage transducer time constant
        VTmin,    ///< \f$V_{\mathrm{Tmin}}\f$ SCL minimum voltage reference value OEL
        VTreset,  ///< \f$V_{\mathrm{Treset}}\f$ SCL voltage reset value b, d OEL
        IQminOEL, ///< \f$I_{\mathrm{QminOEL}}\f$ SCL minimum reactive current reference value OEL
        IQmaxUEL, ///< \f$I_{\mathrm{QmaxUEL}}\f$ SCL maximum reactive current reference value UEL
        KPref     ///< \f$K_{\mathrm{Pref}}\f$ SCL reference scaling factor based on active current
      };

      /// Bus keys for the SCL2C model; deferred until port integration.
      enum class Scl2cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the SCL2C model; deferred until port integration.
      enum class Scl2cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the SCL2C model; deferred until port integration.
      enum class Scl2cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the SCL2C model; deferred until implementation.
      enum class Scl2cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for SCL2C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Scl2cData : public ComponentData<real_type,
                                              index_type,
                                              Scl2cParameters,
                                              Scl2cBuses,
                                              Scl2cSignalInputs,
                                              Scl2cSignalOutputs,
                                              Scl2cMonitorableVariables>
      {
        Scl2cData() = default;

        using Parameters           = Scl2cParameters;
        using Buses                = Scl2cBuses;
        using SignalInputs         = Scl2cSignalInputs;
        using SignalOutputs        = Scl2cSignalOutputs;
        using MonitorableVariables = Scl2cMonitorableVariables;
      };
    } // namespace Limiter
  } // namespace PhasorDynamics
} // namespace GridKit
