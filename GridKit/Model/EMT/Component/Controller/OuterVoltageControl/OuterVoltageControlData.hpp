/**
 * @file OuterVoltageControlData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for the EMT outer-loop voltage controller.
 */
#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      enum class OuterVoltageControlParameters
      {
        C,   ///< \f$C\f$ Filter capacitance [F]
        Kp,  ///< \f$K_P\f$ Proportional gain [S]
        Ki,  ///< \f$K_I\f$ Integral gain [S/s]
        Kaw, ///< \f$K_{\mathrm{aw}}\f$ Tracking anti-windup gain [1/s]
      };

      enum class OuterVoltageControlInputs : size_t
      {
        vrefd, ///< \f$v_d^{\mathrm{ref}}\f$ Voltage reference [V]
        vrefq, ///< \f$v_q^{\mathrm{ref}}\f$ Voltage reference [V]
        vd,    ///< \f$v_d\f$ Capacitor voltage [V]
        vq,    ///< \f$v_q\f$ Capacitor voltage [V]
        igd,   ///< \f$i_{g,d}\f$ Grid-side current [A]
        igq,   ///< \f$i_{g,q}\f$ Grid-side current [A]
        omega, ///< \f$\omega\f$ Electrical angular frequency [rad/s]
        ilimd, ///< \f$i_d^{\mathrm{lim}}\f$ Limited current command [A]
        ilimq, ///< \f$i_q^{\mathrm{lim}}\f$ Limited current command [A]
        SIZE,
      };

      enum class OuterVoltageControlOutputs : size_t
      {
        icmdd, ///< \f$i_d^{\mathrm{cmd}}\f$ Current command [A]
        icmdq, ///< \f$i_q^{\mathrm{cmd}}\f$ Current command [A]
        SIZE,
      };

      enum class OuterVoltageControlMonitorableVariables
      {
        etad,  ///< \f$\eta_d\f$ Integral contribution [A]
        etaq,  ///< \f$\eta_q\f$ Integral contribution [A]
        icmdd, ///< \f$i_d^{\mathrm{cmd}}\f$ Current command [A]
        icmdq, ///< \f$i_q^{\mathrm{cmd}}\f$ Current command [A]
      };

      template <typename real_type, typename index_type>
      struct OuterVoltageControlData : public ComponentData<real_type,
                                                            index_type,
                                                            OuterVoltageControlParameters,
                                                            OuterVoltageControlInputs,
                                                            OuterVoltageControlOutputs,
                                                            OuterVoltageControlMonitorableVariables>
      {
        OuterVoltageControlData() = default;

        using Parameters           = OuterVoltageControlParameters;
        using Inputs               = OuterVoltageControlInputs;
        using Outputs              = OuterVoltageControlOutputs;
        using MonitorableVariables = OuterVoltageControlMonitorableVariables;
      };
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
