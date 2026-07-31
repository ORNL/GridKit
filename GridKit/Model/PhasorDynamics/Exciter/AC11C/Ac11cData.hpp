/**
 * @file Ac11cData.hpp
 * @brief Modeling data for the AC11C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Parameter keys for the AC11C model.
      enum class Ac11cParameters
      {
        RC,     ///< \f$R_{\mathrm{C}}\f$ Resistive component of load compensation
        XC,     ///< \f$X_{\mathrm{C}}\f$ Reactance component of load compensation
        TR,     ///< \f$T_{\mathrm{R}}\f$ Regulator input filter time constant
        KPA,    ///< \f$K_{\mathrm{PA}}\f$ Voltage regulator proportional gain
        TIA,    ///< \f$T_{\mathrm{IA}}\f$ Voltage regulator integral time constant
        KPU,    ///< \f$K_{\mathrm{PU}}\f$ UEL regulator proportional gain
        TIU,    ///< \f$T_{\mathrm{IU}}\f$ UEL regulator integral time constant
        KB,     ///< \f$K_{\mathrm{B}}\f$ Voltage and UEL regulator derivative gain
        TB,     ///< \f$T_{\mathrm{B}}\f$ Time constant for derivative element
        KPO,    ///< \f$K_{\mathrm{PO}}\f$ OEL regulator proportional gain
        TIO,    ///< \f$T_{\mathrm{IO}}\f$ OEL regulator integral time constant
        VRSmax, ///< \f$V_{\mathrm{RSmax}}\f$ Maximum PSS regulator output
        VRSmin, ///< \f$V_{\mathrm{RSmin}}\f$ Minimum PSS regulator output
        VRmax,  ///< \f$V_{\mathrm{Rmax}}\f$ Maximum regulator output
        VRmin,  ///< \f$V_{\mathrm{Rmin}}\f$ Minimum regulator output
        VAmax,  ///< \f$V_{\mathrm{Amax}}\f$ Maximum exciter output
        VAmin,  ///< \f$V_{\mathrm{Amin}}\f$ Minimum exciter output
        TE,     ///< \f$T_{\mathrm{E}}\f$ Exciter field time constant
        KC,     ///< \f$K_{\mathrm{C}}\f$ Rectifier loading factor proportional to commutating reactance
        KD,     ///< \f$K_{\mathrm{D}}\f$ Demagnetizing factor, function of exciter alternator reactances
        KE,     ///< \f$K_{\mathrm{E}}\f$ Exciter field proportional constant
        VFEmax, ///< \f$V_{\mathrm{FEmax}}\f$ Maximum exciter field current
        E1,     ///< \f$E_{1}\f$ Exciter output voltage for saturation factor S (E ) E 1
        Se1,    ///< \f$S_{\mathrm{E}}(E_{1})\f$ Exciter saturation factor at exciter output voltage E 1
        E2,     ///< \f$E_{2}\f$ Exciter output voltage for saturation factor S (E ) E 2
        Se2,    ///< \f$S_{\mathrm{E}}(E_{2})\f$ Exciter saturation factor at exciter output voltage E 2
        SW1,    ///< \f$\mathrm{SW}_{1}\f$ Power source selector
        KP,     ///< \f$K_{\mathrm{P}}\f$ Potential circuit (voltage) gain coefficient
        KI,     ///< \f$K_{\mathrm{I}}\f$ Potential circuit (current) gain coefficient
        XL,     ///< \f$X_{\mathrm{L}}\f$ Reactance associated with potential source
        ThetaP, ///< \f$\theta_{\mathrm{P}}\f$ Potential circuit phase angle (degrees)
        KC1,    ///< \f$K_{\mathrm{C1}}\f$ Rectifier loading factor proportional to commutating reactance
        VB1max, ///< \f$V_{\mathrm{B1max}}\f$ Maximum available exciter voltage
        KI2,    ///< \f$K_{\mathrm{I2}}\f$ Additive potential circuit (current) gain coefficient
        KC2,    ///< \f$K_{\mathrm{C2}}\f$ Rectifier loading factor proportional to commutating reactance
        VB2max, ///< \f$V_{\mathrm{B2max}}\f$ Maximum available additive exciter voltage
        KBOOST, ///< \f$K_{\mathrm{BOOST}}\f$ Additive independent source
        VBOOST  ///< \f$V_{\mathrm{BOOST}}\f$ Reference value for applying additive (boost) circuit [see note 2]
      };

      /// Bus keys for the AC11C model; deferred until port integration.
      enum class Ac11cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the AC11C model; deferred until port integration.
      enum class Ac11cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the AC11C model; deferred until port integration.
      enum class Ac11cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the AC11C model; deferred until implementation.
      enum class Ac11cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for AC11C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Ac11cData : public ComponentData<real_type,
                                              index_type,
                                              Ac11cParameters,
                                              Ac11cBuses,
                                              Ac11cSignalInputs,
                                              Ac11cSignalOutputs,
                                              Ac11cMonitorableVariables>
      {
        Ac11cData() = default;

        using Parameters           = Ac11cParameters;
        using Buses                = Ac11cBuses;
        using SignalInputs         = Ac11cSignalInputs;
        using SignalOutputs        = Ac11cSignalOutputs;
        using MonitorableVariables = Ac11cMonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
