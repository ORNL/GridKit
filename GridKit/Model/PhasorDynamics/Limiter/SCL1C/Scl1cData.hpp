/**
 * @file Scl1cData.hpp
 * @brief Modeling data for the SCL1C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Limiter
    {
      /// Parameter keys for the SCL1C model.
      enum class Scl1cParameters
      {
        ISCLlim, ///< \f$I_{\mathrm{SCLlim}}\f$ SCL terminal current pick up level
        TIT,     ///< \f$T_{\mathrm{IT}}\f$ Terminal current transducer equivalent time constant
        K,       ///< \f$K\f$ SCL timing characteristic factor
        TQSCL,   ///< \f$T_{\mathrm{QSCL}}\f$ Reactive current transducer equivalent time constant
        IQmin,   ///< \f$I_{\mathrm{Qmin}}\f$ Dead-band for reactive current
        VSCLdb,  ///< \f$V_{\mathrm{SCLdb}}\f$ Dead-band for reactive power or power factor
        TINV,    ///< \f$T_{\mathrm{INV}}\f$ Inverse time delay after pickup
        TDSCL,   ///< \f$T_{\mathrm{DSCL}}\f$ Fixed-time delay after pickup
        SW1,     ///< \f$\mathrm{SW}_{1}\f$ Reactive current/reactive power selector
        SW2,     ///< \f$\mathrm{SW}_{2}\f$ Fixed-time or inverse time selector
        KPoex,   ///< \f$K_{\mathrm{Poex}}\f$ SCL proportional gain (overexcited range)
        KIoex,   ///< \f$K_{\mathrm{Ioex}}\f$ SCL integral gain (overexcited range)
        KPuex,   ///< \f$K_{\mathrm{Puex}}\f$ SCL proportional gain (underexcited range)
        KIuex,   ///< \f$K_{\mathrm{Iuex}}\f$ SCL integral gain (underexcited range)
        VSCLmax, ///< \f$V_{\mathrm{SCLmax}}\f$ SCL upper integrator limit
        VSCLmin  ///< \f$V_{\mathrm{SCLmin}}\f$ SCL lower integrator limit
      };

      /// Bus keys for the SCL1C model; deferred until port integration.
      enum class Scl1cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the SCL1C model; deferred until port integration.
      enum class Scl1cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the SCL1C model; deferred until port integration.
      enum class Scl1cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the SCL1C model; deferred until implementation.
      enum class Scl1cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for SCL1C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Scl1cData : public ComponentData<real_type,
                                              index_type,
                                              Scl1cParameters,
                                              Scl1cBuses,
                                              Scl1cSignalInputs,
                                              Scl1cSignalOutputs,
                                              Scl1cMonitorableVariables>
      {
        Scl1cData() = default;

        using Parameters           = Scl1cParameters;
        using Buses                = Scl1cBuses;
        using SignalInputs         = Scl1cSignalInputs;
        using SignalOutputs        = Scl1cSignalOutputs;
        using MonitorableVariables = Scl1cMonitorableVariables;
      };
    } // namespace Limiter
  } // namespace PhasorDynamics
} // namespace GridKit
