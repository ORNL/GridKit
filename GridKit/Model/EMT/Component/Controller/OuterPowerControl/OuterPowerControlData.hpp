/**
 * @file OuterPowerControlData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for the EMT outer power controller.
 */
#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      enum class OuterPowerControlParameters
      {
        I,   ///< Nominal phase RMS current for absolute tolerances [A]
        V,   ///< \f$V\f$ Rated line-to-line RMS voltage [V]
        Kp,  ///< \f$K_P\f$ Proportional gain [-]
        Ki,  ///< \f$K_I\f$ Integral gain [1/s]
        Kaw, ///< \f$K_{\mathrm{aw}}\f$ Tracking anti-windup gain [1/s]
      };

      enum class OuterPowerControlInputs : size_t
      {
        vd,    ///< \f$v_d\f$ Terminal voltage [V]
        vq,    ///< \f$v_q\f$ Terminal voltage [V]
        id,    ///< \f$i_d\f$ Measured current [A]
        iq,    ///< \f$i_q\f$ Measured current [A]
        ilimd, ///< \f$i_d^{\mathrm{lim}}\f$ Limited current command [A]
        ilimq, ///< \f$i_q^{\mathrm{lim}}\f$ Limited current command [A]
        Pref,  ///< \f$P^{\mathrm{ref}}\f$ Active-power reference [W]
        Qref,  ///< \f$Q^{\mathrm{ref}}\f$ Reactive-power reference [var]
        SIZE,
      };

      enum class OuterPowerControlOutputs : size_t
      {
        icmdd, ///< \f$i_d^{\mathrm{cmd}}\f$ Current command [A]
        icmdq, ///< \f$i_q^{\mathrm{cmd}}\f$ Current command [A]
        SIZE,
      };

      enum class OuterPowerControlMonitorableVariables
      {
        etad,  ///< \f$\eta_d\f$ Integral contribution [A]
        etaq,  ///< \f$\eta_q\f$ Integral contribution [A]
        icmdd, ///< \f$i_d^{\mathrm{cmd}}\f$ Current command [A]
        icmdq, ///< \f$i_q^{\mathrm{cmd}}\f$ Current command [A]
      };

      template <typename real_type, typename index_type>
      struct OuterPowerControlData : public ComponentData<real_type,
                                                          index_type,
                                                          OuterPowerControlParameters,
                                                          OuterPowerControlInputs,
                                                          OuterPowerControlOutputs,
                                                          OuterPowerControlMonitorableVariables>
      {
        OuterPowerControlData() = default;

        using Parameters           = OuterPowerControlParameters;
        using Inputs               = OuterPowerControlInputs;
        using Outputs              = OuterPowerControlOutputs;
        using MonitorableVariables = OuterPowerControlMonitorableVariables;
      };
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
