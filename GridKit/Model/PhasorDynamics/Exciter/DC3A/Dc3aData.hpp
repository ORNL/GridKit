/**
 * @file Dc3aData.hpp
 * @brief Modeling data for the DC3A model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Parameter keys for the DC3A model.
      enum class Dc3aParameters
      {
        RC,    ///< \f$R_{\mathrm{C}}\f$ Resistive component of load compensation
        XC,    ///< \f$X_{\mathrm{C}}\f$ Reactance component of load compensation
        TR,    ///< \f$T_{\mathrm{R}}\f$ Regulator input filter time constant
        KE,    ///< \f$K_{\mathrm{E}}\f$ Exciter field proportional constant
        TE,    ///< \f$T_{\mathrm{E}}\f$ Exciter field time constant
        KV,    ///< \f$K_{\mathrm{V}}\f$ Fast raise/lower contact setting
        VRmax, ///< \f$V_{\mathrm{Rmax}}\f$ Maximum controller output
        VRmin, ///< \f$V_{\mathrm{Rmin}}\f$ Minimum controller output
        TRH,   ///< \f$T_{\mathrm{RH}}\f$ Rheostat travel time
        E1,    ///< \f$E_{1}\f$ Exciter output voltage for saturation factor S (E ) E 1
        SE1,   ///< \f$S_{\mathrm{E1}}\f$ Exciter saturation factor at exciter output voltage E 1
        E2,    ///< \f$E_{2}\f$ Exciter output voltage for saturation factor S (E ) E 2
        SE2    ///< \f$S_{\mathrm{E2}}\f$ Exciter saturation factor at exciter output voltage E 2
      };

      /// Bus keys for the DC3A model; deferred until port integration.
      enum class Dc3aBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the DC3A model; deferred until port integration.
      enum class Dc3aSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the DC3A model; deferred until port integration.
      enum class Dc3aSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the DC3A model; deferred until implementation.
      enum class Dc3aMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for DC3A.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Dc3aData : public ComponentData<real_type,
                                             index_type,
                                             Dc3aParameters,
                                             Dc3aBuses,
                                             Dc3aSignalInputs,
                                             Dc3aSignalOutputs,
                                             Dc3aMonitorableVariables>
      {
        Dc3aData() = default;

        using Parameters           = Dc3aParameters;
        using Buses                = Dc3aBuses;
        using SignalInputs         = Dc3aSignalInputs;
        using SignalOutputs        = Dc3aSignalOutputs;
        using MonitorableVariables = Dc3aMonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
