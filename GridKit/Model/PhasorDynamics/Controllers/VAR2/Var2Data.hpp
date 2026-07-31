/**
 * @file Var2Data.hpp
 * @brief Modeling data for the VAR2 model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Controller
    {
      /// Parameter keys for the VAR2 model.
      enum class Var2Parameters
      {
        QREF,   ///< \f$Q_{\mathrm{REF}}\f$ Var controller reference setpoint
        VITmin, ///< \f$V_{\mathrm{ITmin}}\f$ Var controller minimum terminal current limit
        VVTmin, ///< \f$V_{\mathrm{VTmin}}\f$ Var controller minimum terminal voltage limit
        VVTmax, ///< \f$V_{\mathrm{VTmax}}\f$ Var controller maximum terminal voltage limit
        KPvar,  ///< \f$K_{\mathrm{Pvar}}\f$ Var controller proportional gain
        KIvar,  ///< \f$K_{\mathrm{Ivar}}\f$ Var controller integral gain
        VVARLMT ///< \f$V_{\mathrm{VARLMT}}\f$ Var controller output limit
      };

      /// Bus keys for the VAR2 model; deferred until port integration.
      enum class Var2Buses : size_t
      {
        SIZE
      };

      /// Signal input keys for the VAR2 model; deferred until port integration.
      enum class Var2SignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the VAR2 model; deferred until port integration.
      enum class Var2SignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the VAR2 model; deferred until implementation.
      enum class Var2MonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for VAR2.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Var2Data : public ComponentData<real_type,
                                             index_type,
                                             Var2Parameters,
                                             Var2Buses,
                                             Var2SignalInputs,
                                             Var2SignalOutputs,
                                             Var2MonitorableVariables>
      {
        Var2Data() = default;

        using Parameters           = Var2Parameters;
        using Buses                = Var2Buses;
        using SignalInputs         = Var2SignalInputs;
        using SignalOutputs        = Var2SignalOutputs;
        using MonitorableVariables = Var2MonitorableVariables;
      };
    } // namespace Controller
  } // namespace PhasorDynamics
} // namespace GridKit
