/**
 * @file GasBranchData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for gas branches (pipes).
 */
#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Gas
    {
      /// Parameters for a gas branch
      enum class GasBranchParameters : size_t
      {
        L,      ///< \f$L\f$ Pipe length [m]
        D,      ///< \f$D\f$ Pipe diameter [m]
        lambda, ///< \f$\lambda\f$ Darcy friction factor [-]
        a,      ///< \f$a\f$ Isothermal speed of sound [m/s]
        b1,     ///< \f$b_1\f$ First density coefficient [-]
        b2,     ///< \f$b_2\f$ Second density coefficient [1/Pa]
        dxmax,  ///< \f$\Delta x^{\max}\f$ Maximum segment length [m]
      };

      /// Buses for a gas branch
      enum class GasBranchBuses : size_t
      {
        bus1, ///< \f$p_0\f$ Required Known bus-1 terminal pressure [Pa]
        bus2, ///< \f$p_N\f$ Required Known bus-2 terminal pressure [Pa]
      };

      /// Signal inputs supported for a gas branch
      enum class GasBranchSignalInputs : size_t
      {
      };

      /// Signal outputs supported for a gas branch
      enum class GasBranchSignalOutputs : size_t
      {
      };

      /// Variables able to be monitored for a gas branch
      enum class GasBranchMonitorableVariables : size_t
      {
        q1,       ///< \f$\Delta q^{\mathrm{bus}}_1\f$ Bus-1 terminal mass flow, entering bus 1 [kg/s]
        q2,       ///< \f$\Delta q^{\mathrm{bus}}_2\f$ Bus-2 terminal mass flow, entering bus 2 [kg/s]
        linepack, ///< \f$m\f$ Gas mass in the pipe [kg]
      };

      /**
       * @brief Contains modeling data for a GasBranch
       *
       * @tparam real_type  Real parameter data type
       * @tparam index_type Integer parameter data type
       */
      template <typename real_type, typename index_type>
      using GasBranchData =
          ComponentData<real_type,
                        index_type,
                        GasBranchParameters,
                        GasBranchBuses,
                        GasBranchSignalInputs,
                        GasBranchSignalOutputs,
                        GasBranchMonitorableVariables>;
    } // namespace Gas
  } // namespace PhasorDynamics
} // namespace GridKit
