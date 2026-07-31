/**
 * @file Oel2cData.hpp
 * @brief Modeling data for the OEL2C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Limiter
    {
      /// Parameter keys for the OEL2C model.
      enum class Oel2cParameters
      {
        TC1oel,   ///< \f$T_{\mathrm{C1oel}}\f$ OEL regulator denominator (lag) time constant 1
        TB1oel,   ///< \f$T_{\mathrm{B1oel}}\f$ OEL regulator numerator (lead) time constant 1
        TC2oel,   ///< \f$T_{\mathrm{C2oel}}\f$ OEL regulator denominator (lag) time constant 2
        TB2oel,   ///< \f$T_{\mathrm{B2oel}}\f$ OEL regulator numerator (lead) time constant 2
        KPoel,    ///< \f$K_{\mathrm{Poel}}\f$ OEL PID regulator proportional gain
        KIoel,    ///< \f$K_{\mathrm{Ioel}}\f$ OEL PID regulator integral gain
        KDoel,    ///< \f$K_{\mathrm{Doel}}\f$ OEL PID regulator differential gain
        TDoel,    ///< \f$T_{\mathrm{Doel}}\f$ OEL PID regulator differential time constant
        VOELmax3, ///< \f$V_{\mathrm{OELmax3}}\f$ Maximum OEL PID output limit
        VOELmin3, ///< \f$V_{\mathrm{OELmin3}}\f$ Minimum OEL PID output limit
        VOELmax2, ///< \f$V_{\mathrm{OELmax2}}\f$ Maximum OEL lead-lag 1 output limit
        VOELmin2, ///< \f$V_{\mathrm{OELmin2}}\f$ Minimum OEL lead-lag 1 output limit
        VOELmax1, ///< \f$V_{\mathrm{OELmax1}}\f$ Maximum OEL output limit
        VOELmin1, ///< \f$V_{\mathrm{OELmin1}}\f$ Minimum OEL output limit
        Ireset,   ///< \f$I_{\mathrm{reset}}\f$ OEL reset-reference, if OEL is inactive
        Ten,      ///< \f$T_{\mathrm{en}}\f$ OEL activation delay time
        Toff,     ///< \f$T_{\mathrm{off}}\f$ OEL reset delay time
        ITHoff,   ///< \f$I_{\mathrm{THoff}}\f$ OEL reset threshold value
        KSCALE,   ///< \f$K_{\mathrm{SCALE}}\f$ OEL input signal scaling factord
        TRoel,    ///< \f$T_{\mathrm{Roel}}\f$ OEL input signal filter time constant
        Kact,     ///< \f$K_{\mathrm{act}}\f$ OEL actual value scaling factor
        ITFpu,    ///< \f$I_{\mathrm{TFpu}}\f$ OEL reference for inverse time calculations
        Iinst,    ///< \f$I_{\mathrm{inst}}\f$ OEL instantaneous field current limit
        Ilim,     ///< \f$I_{\mathrm{lim}}\f$ OEL thermal field current limit
        TAoel,    ///< \f$T_{\mathrm{Aoel}}\f$ OEL reference filter time constant
        c1,       ///< \f$c_{1}\f$ OEL exponent for calculation of I ERRinv1
        K1,       ///< \f$K_{1}\f$ OEL gain for calculation of I ERRinv1
        c2,       ///< \f$c_{2}\f$ OEL exponent for calculation of I b ERRinv2
        K2,       ///< \f$K_{2}\f$ OEL gain for calculation of I b ERRinv2
        VINVmax,  ///< \f$V_{\mathrm{INVmax}}\f$ OEL maximum inverse time output
        VINVmin,  ///< \f$V_{\mathrm{INVmin}}\f$ OEL minimum inverse time output
        Fixedru,  ///< \f$\mathrm{Fixed}_{\mathrm{ru}}\f$ OEL fixed delay time output
        Fixedrd,  ///< \f$\mathrm{Fixed}_{\mathrm{rd}}\f$ OEL fixed cooling-down time output
        TFCL,     ///< \f$T_{\mathrm{FCL}}\f$ OEL timer reference
        Tmax,     ///< \f$T_{\mathrm{max}}\f$ OEL timer maximum level
        Tmin,     ///< \f$T_{\mathrm{min}}\f$ OEL timer minimum level
        KFB,      ///< \f$K_{\mathrm{FB}}\f$ OEL timer feedback gain
        Krd,      ///< \f$K_{\mathrm{rd}}\f$ OEL reference ramp-down rate
        Kru,      ///< \f$K_{\mathrm{ru}}\f$ OEL reference ramp-up rate
        KZRU,     ///< \f$K_{\mathrm{ZRU}}\f$ OEL thermal reference release threshold
        IFDrated  ///< \f$I_{\mathrm{FDrated}}\f$ Rated field current
      };

      /// Bus keys for the OEL2C model; deferred until port integration.
      enum class Oel2cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the OEL2C model; deferred until port integration.
      enum class Oel2cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the OEL2C model; deferred until port integration.
      enum class Oel2cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the OEL2C model; deferred until implementation.
      enum class Oel2cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for OEL2C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Oel2cData : public ComponentData<real_type,
                                              index_type,
                                              Oel2cParameters,
                                              Oel2cBuses,
                                              Oel2cSignalInputs,
                                              Oel2cSignalOutputs,
                                              Oel2cMonitorableVariables>
      {
        Oel2cData() = default;

        using Parameters           = Oel2cParameters;
        using Buses                = Oel2cBuses;
        using SignalInputs         = Oel2cSignalInputs;
        using SignalOutputs        = Oel2cSignalOutputs;
        using MonitorableVariables = Oel2cMonitorableVariables;
      };
    } // namespace Limiter
  } // namespace PhasorDynamics
} // namespace GridKit
