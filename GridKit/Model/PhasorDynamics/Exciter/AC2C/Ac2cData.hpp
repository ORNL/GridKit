/**
 * @file Ac2cData.hpp
 * @brief Modeling data for the AC2C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Parameter keys for the AC2C model.
      enum class Ac2cParameters
      {
        RC,     ///< \f$R_{\mathrm{C}}\f$ Resistive component of load compensation
        XC,     ///< \f$X_{\mathrm{C}}\f$ Reactance component of load compensation
        TR,     ///< \f$T_{\mathrm{R}}\f$ Regulator input filter time constant
        KA,     ///< \f$K_{\mathrm{A}}\f$ Regulator output gain
        TA,     ///< \f$T_{\mathrm{A}}\f$ Regulator output time constant
        TB,     ///< \f$T_{\mathrm{B}}\f$ Regulator denominator (lag) time constant
        TC,     ///< \f$T_{\mathrm{C}}\f$ Regulator numerator (lead) time constant
        KB,     ///< \f$K_{\mathrm{B}}\f$ Second stage regulator gain
        KH,     ///< \f$K_{\mathrm{H}}\f$ Exciter field current regulator feedback gain
        KF,     ///< \f$K_{\mathrm{F}}\f$ Rate feedback excitation system stabilizer gain
        TF,     ///< \f$T_{\mathrm{F}}\f$ Rate feedback time constant
        KE,     ///< \f$K_{\mathrm{E}}\f$ Exciter field proportional constant
        TE,     ///< \f$T_{\mathrm{E}}\f$ Exciter field time constant
        KD,     ///< \f$K_{\mathrm{D}}\f$ Demagnetizing factor, function of exciter alternator reactances
        KC,     ///< \f$K_{\mathrm{C}}\f$ Rectifier loading factor proportional to commutating reactance
        VAmax,  ///< \f$V_{\mathrm{Amax}}\f$ Maximum regulator output
        VAmin,  ///< \f$V_{\mathrm{Amin}}\f$ Minimum regulator output
        EFEmax, ///< \f$E_{\mathrm{FEmax}}\f$ Maximum exciter field voltage
        EFEmin, ///< \f$E_{\mathrm{FEmin}}\f$ Minimum exciter field voltage
        VFEmax, ///< \f$V_{\mathrm{FEmax}}\f$ Maximum exciter field current limit reference
        VEmin,  ///< \f$V_{\mathrm{Emin}}\f$ Minimum exciter voltage output
        E1,     ///< \f$E_{1}\f$ Exciter output voltage for saturation factor S (E ) E 1
        Se1,    ///< \f$S_{\mathrm{E}}(E_{1})\f$ Exciter saturation factor at exciter output voltage E 1
        E2,     ///< \f$E_{2}\f$ Exciter output voltage for saturation factor S (E ) E 2
        Se2     ///< \f$S_{\mathrm{E}}(E_{2})\f$ Exciter saturation factor at exciter output voltage E 2
      };

      /// Bus keys for the AC2C model; deferred until port integration.
      enum class Ac2cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the AC2C model; deferred until port integration.
      enum class Ac2cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the AC2C model; deferred until port integration.
      enum class Ac2cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the AC2C model; deferred until implementation.
      enum class Ac2cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for AC2C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Ac2cData : public ComponentData<real_type,
                                             index_type,
                                             Ac2cParameters,
                                             Ac2cBuses,
                                             Ac2cSignalInputs,
                                             Ac2cSignalOutputs,
                                             Ac2cMonitorableVariables>
      {
        Ac2cData() = default;

        using Parameters           = Ac2cParameters;
        using Buses                = Ac2cBuses;
        using SignalInputs         = Ac2cSignalInputs;
        using SignalOutputs        = Ac2cSignalOutputs;
        using MonitorableVariables = Ac2cMonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
