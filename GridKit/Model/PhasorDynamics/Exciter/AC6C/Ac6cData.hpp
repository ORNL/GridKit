/**
 * @file Ac6cData.hpp
 * @brief Modeling data for the AC6C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Parameter keys for the AC6C model.
      enum class Ac6cParameters
      {
        RC,     ///< \f$R_{\mathrm{C}}\f$ Resistive component of load compensation
        XC,     ///< \f$X_{\mathrm{C}}\f$ Reactance component of load compensation
        TR,     ///< \f$T_{\mathrm{R}}\f$ Regulator input filter time constant
        KA,     ///< \f$K_{\mathrm{A}}\f$ Regulator output gain
        TA,     ///< \f$T_{\mathrm{A}}\f$ Regulator output time constant
        TC,     ///< \f$T_{\mathrm{C}}\f$ Regulator numerator (lead) time constant
        TB,     ///< \f$T_{\mathrm{B}}\f$ Regulator denominator (lag) time constant
        TK,     ///< \f$T_{\mathrm{K}}\f$ Regulator numerator (lead) time constant
        KH,     ///< \f$K_{\mathrm{H}}\f$ Exciter field current limiter gain
        TE,     ///< \f$T_{\mathrm{E}}\f$ Exciter field time constant
        TH,     ///< \f$T_{\mathrm{H}}\f$ Exciter field current limiter denominator (lag) time constant
        TJ,     ///< \f$T_{\mathrm{J}}\f$ Exciter field current limiter numerator (lead) time constant
        KC,     ///< \f$K_{\mathrm{C}}\f$ Rectifier loading factor proportional to commutating reactance
        KD,     ///< \f$K_{\mathrm{D}}\f$ Demagnetizing factor, function of exciter alternator reactances
        KE,     ///< \f$K_{\mathrm{E}}\f$ Exciter field proportional constant
        VAmax,  ///< \f$V_{\mathrm{Amax}}\f$ Maximum voltage regulator outputs
        VAmin,  ///< \f$V_{\mathrm{Amin}}\f$ Minimum voltage regulator outputs
        EFEmax, ///< \f$E_{\mathrm{FEmax}}\f$ Maximum exciter field voltage
        EFEmin, ///< \f$E_{\mathrm{FEmin}}\f$ Minimum exciter field voltage
        VHmax,  ///< \f$V_{\mathrm{Hmax}}\f$ Exciter field current limiter maximum output
        VFELIM, ///< \f$V_{\mathrm{FELIM}}\f$ Exciter field current limiter reference
        VFEmax, ///< \f$V_{\mathrm{FEmax}}\f$ Maximum exciter field current
        E1,     ///< \f$E_{1}\f$ Exciter output voltage for saturation factor S (E ) E 1
        Se1,    ///< \f$S_{\mathrm{E}}(E_{1})\f$ Exciter saturation factor at exciter output voltage E 1
        E2,     ///< \f$E_{2}\f$ Exciter output voltage for saturation factor S (E ) E 2
        Se2     ///< \f$S_{\mathrm{E}}(E_{2})\f$ Exciter saturation factor at exciter output voltage E 2
      };

      /// Bus keys for the AC6C model; deferred until port integration.
      enum class Ac6cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the AC6C model; deferred until port integration.
      enum class Ac6cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the AC6C model; deferred until port integration.
      enum class Ac6cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the AC6C model; deferred until implementation.
      enum class Ac6cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for AC6C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Ac6cData : public ComponentData<real_type,
                                             index_type,
                                             Ac6cParameters,
                                             Ac6cBuses,
                                             Ac6cSignalInputs,
                                             Ac6cSignalOutputs,
                                             Ac6cMonitorableVariables>
      {
        Ac6cData() = default;

        using Parameters           = Ac6cParameters;
        using Buses                = Ac6cBuses;
        using SignalInputs         = Ac6cSignalInputs;
        using SignalOutputs        = Ac6cSignalOutputs;
        using MonitorableVariables = Ac6cMonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
