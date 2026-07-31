/**
 * @file Pf2Data.hpp
 * @brief Modeling data for the PF2 model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Controller
    {
      /// Parameter keys for the PF2 model.
      enum class Pf2Parameters
      {
        PFREFnorm, ///< \f$\mathrm{PF}_{\mathrm{REFnorm}}\f$ Power factor controller normalized reference setpoint
        VITmin,    ///< \f$V_{\mathrm{ITmin}}\f$ Power factor controller minimum terminal current limit
        VVTmin,    ///< \f$V_{\mathrm{VTmin}}\f$ Power factor controller minimum terminal voltage limit
        VVTmax,    ///< \f$V_{\mathrm{VTmax}}\f$ Power factor controller maximum terminal voltage limit
        KPpf,      ///< \f$K_{\mathrm{Ppf}}\f$ Power factor controller proportional gain
        KIpf,      ///< \f$K_{\mathrm{Ipf}}\f$ Power factor controller integral gain
        VPFLMT     ///< \f$V_{\mathrm{PFLMT}}\f$ Power factor controller output limit
      };

      /// Bus keys for the PF2 model; deferred until port integration.
      enum class Pf2Buses : size_t
      {
        SIZE
      };

      /// Signal input keys for the PF2 model; deferred until port integration.
      enum class Pf2SignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the PF2 model; deferred until port integration.
      enum class Pf2SignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the PF2 model; deferred until implementation.
      enum class Pf2MonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for PF2.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Pf2Data : public ComponentData<real_type,
                                            index_type,
                                            Pf2Parameters,
                                            Pf2Buses,
                                            Pf2SignalInputs,
                                            Pf2SignalOutputs,
                                            Pf2MonitorableVariables>
      {
        Pf2Data() = default;

        using Parameters           = Pf2Parameters;
        using Buses                = Pf2Buses;
        using SignalInputs         = Pf2SignalInputs;
        using SignalOutputs        = Pf2SignalOutputs;
        using MonitorableVariables = Pf2MonitorableVariables;
      };
    } // namespace Controller
  } // namespace PhasorDynamics
} // namespace GridKit
