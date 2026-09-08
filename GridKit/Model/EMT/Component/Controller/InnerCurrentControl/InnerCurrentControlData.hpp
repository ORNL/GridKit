/**
 * @file InnerCurrentControlData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for the EMT inner-loop current controller.
 */
#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      enum class InnerCurrentControlParameters
      {
        L,    ///< \f$L\f$ Filter inductance [H]
        Kp,   ///< \f$K_P\f$ Proportional gain [ohm]
        Ki,   ///< \f$K_I\f$ Integral gain [ohm/s]
        Kaw,  ///< \f$K_{\mathrm{aw}}\f$ Tracking anti-windup gain [1/s]
        Imax, ///< \f$I^{\max}\f$ Current-command norm limit [A]
        Mmax, ///< \f$M^{\max}\f$ Modulation limit [-]
      };

      enum class InnerCurrentControlInputs : size_t
      {
        vd,    ///< \f$v_d\f$ Capacitor voltage [V]
        vq,    ///< \f$v_q\f$ Capacitor voltage [V]
        id,    ///< \f$i_d\f$ Inverter-side filter current [A]
        iq,    ///< \f$i_q\f$ Inverter-side filter current [A]
        icmdd, ///< \f$i_d^{\mathrm{cmd}}\f$ Current command [A]
        icmdq, ///< \f$i_q^{\mathrm{cmd}}\f$ Current command [A]
        omega, ///< \f$\omega\f$ Electrical angular frequency [rad/s]
        vdc,   ///< \f$v_{\mathrm{dc}}\f$ DC-link voltage [V]
        SIZE,
      };

      enum class InnerCurrentControlOutputs : size_t
      {
        ilimd, ///< \f$i_d^{\mathrm{lim}}\f$ Limited current command [A]
        ilimq, ///< \f$i_q^{\mathrm{lim}}\f$ Limited current command [A]
        ud,    ///< \f$u_d\f$ Converter voltage command [V]
        uq,    ///< \f$u_q\f$ Converter voltage command [V]
        SIZE,
      };

      enum class InnerCurrentControlMonitorableVariables
      {
        xid,   ///< \f$\xi_d\f$ Integral contribution [V]
        xiq,   ///< \f$\xi_q\f$ Integral contribution [V]
        ilimd, ///< \f$i_d^{\mathrm{lim}}\f$ Limited current command [A]
        ilimq, ///< \f$i_q^{\mathrm{lim}}\f$ Limited current command [A]
        ud,    ///< \f$u_d\f$ Converter voltage command [V]
        uq,    ///< \f$u_q\f$ Converter voltage command [V]
      };

      template <typename real_type, typename index_type>
      struct InnerCurrentControlData : public ComponentData<real_type,
                                                            index_type,
                                                            InnerCurrentControlParameters,
                                                            InnerCurrentControlInputs,
                                                            InnerCurrentControlOutputs,
                                                            InnerCurrentControlMonitorableVariables>
      {
        InnerCurrentControlData() = default;

        using Parameters           = InnerCurrentControlParameters;
        using Inputs               = InnerCurrentControlInputs;
        using Outputs              = InnerCurrentControlOutputs;
        using MonitorableVariables = InnerCurrentControlMonitorableVariables;
      };
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
