/**
 * @file Oel5cData.hpp
 * @brief Modeling data for the OEL5C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Limiter
    {
      /// Parameter keys for the OEL5C model.
      enum class Oel5cParameters
      {
        IFDpu,    ///< \f$I_{\mathrm{FDpu}}\f$ OEL inverse time integrator pickup level
        IFDlim,   ///< \f$I_{\mathrm{FDlim}}\f$ OEL inverse time limit active level
        VOELmax1, ///< \f$V_{\mathrm{OELmax1}}\f$ OEL inverse time time upper limit
        TOEL,     ///< \f$T_{\mathrm{OEL}}\f$ OEL inverse time integrator time constant
        KIFDT,    ///< \f$K_{\mathrm{IFDT}}\f$ OEL inverse time leak gain
        K,        ///< \f$K\f$ OEL lead-lag gain
        TCoel,    ///< \f$T_{\mathrm{Coel}}\f$ OEL lead time constant
        TBoel,    ///< \f$T_{\mathrm{Boel}}\f$ OEL lag time constant
        IFDpulev, ///< \f$I_{\mathrm{FDpulev}}\f$ OEL activation logic pickup level
        TIFDlev,  ///< \f$T_{\mathrm{IFDlev}}\f$ OEL activation logic timer setpoint
        IFDref1,  ///< \f$I_{\mathrm{FDref1}}\f$ OEL reference 1
        IFDref2,  ///< \f$I_{\mathrm{FDref2}}\f$ OEL reference 2
        KPoel,    ///< \f$K_{\mathrm{Poel}}\f$ OEL proportional gain
        KIoel,    ///< \f$K_{\mathrm{Ioel}}\f$ OEL integral gain
        VOELmax,  ///< \f$V_{\mathrm{OELmax}}\f$ OEL PI control upper limit
        VOELmin,  ///< \f$V_{\mathrm{OELmin}}\f$ OEL PI control lower limit
        KPvfe,    ///< \f$K_{\mathrm{Pvfe}}\f$ Exciter field current regulator proportional gain
        KIvfe,    ///< \f$K_{\mathrm{Ivfe}}\f$ Exciter field current regulator integral gain
        VVFEmax,  ///< \f$V_{\mathrm{VFEmax}}\f$ Exciter field current regulator upper limit
        VVFEmin,  ///< \f$V_{\mathrm{VFEmin}}\f$ Exciter field current regulator lower limit
        KSCALE1,  ///< \f$K_{\mathrm{SCALE1}}\f$ Scale factor for OEL input
        TF2,      ///< \f$T_{\mathrm{F2}}\f$ OEL input or exciter field current transducer time constant
        KSCALE2,  ///< \f$K_{\mathrm{SCALE2}}\f$ Scale factor I /I FEbase FErated
        VFEref,   ///< \f$V_{\mathrm{FEref}}\f$ Exciter field current reference setpoint
        SW1,      ///< \f$\mathrm{SW}_{1}\f$ OEL reference logic switch
        Ibias,    ///< \f$I_{\mathrm{bias}}\f$ OEL reference bias
        K1        ///< \f$K_{1}\f$ Exponent for inverse time function
      };

      /// Bus keys for the OEL5C model; deferred until port integration.
      enum class Oel5cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the OEL5C model; deferred until port integration.
      enum class Oel5cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the OEL5C model; deferred until port integration.
      enum class Oel5cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the OEL5C model; deferred until implementation.
      enum class Oel5cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for OEL5C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Oel5cData : public ComponentData<real_type,
                                              index_type,
                                              Oel5cParameters,
                                              Oel5cBuses,
                                              Oel5cSignalInputs,
                                              Oel5cSignalOutputs,
                                              Oel5cMonitorableVariables>
      {
        Oel5cData() = default;

        using Parameters           = Oel5cParameters;
        using Buses                = Oel5cBuses;
        using SignalInputs         = Oel5cSignalInputs;
        using SignalOutputs        = Oel5cSignalOutputs;
        using MonitorableVariables = Oel5cMonitorableVariables;
      };
    } // namespace Limiter
  } // namespace PhasorDynamics
} // namespace GridKit
