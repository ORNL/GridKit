/**
 * @file GensalData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for a GENSAL generator model.
 *
 */
#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Initial parameters for a Gensal generator model
    enum class GensalParameters
    {
      p0,    ///< \f$P_0\f$ Initial active power
      q0,    ///< \f$Q_0\f$ Initial reactive power
      H,     ///< \f$H\f$ Rotor inertia
      D,     ///< \f$D\f$ Damping coefficient
      Ra,    ///< \f$R_a\f$ Winding resistance
      Tdop,  ///< \f$T'_{d0}\f$ Open circuit direct axis transient time
      Tdopp, ///< \f$T''_{d0}\f$ Open circuit direct axis sub-transient time
      Tqopp, ///< \f$T''_{q0}\f$ Open circuit quadrature axis sub-transient time
      Xd,    ///< \f$X_d\f$ Direct axis synchronous reactance
      Xdp,   ///< \f$X'_d\f$ Direct axis transient reactance
      Xdpp,  ///< \f$X''_d\f$ Direct axis sub-transient reactance
      Xq,    ///< \f$X_q\f$ Quadrature axis synchronous reactance
      Xl,    ///< \f$X_{\ell}\f$ Stator leakage reactance
      S10,   ///< \f$S_{10}\f$ Saturation factor at 1.0 pu flux
      S12,   ///< \f$S_{12}\f$ Saturation factor at 1.2 pu flux
      mva,   ///< \f$S_\mathrm{mach}\f$ MVA base of the gensal model
    };

    /// Buses for a Gensal generator model
    enum class GensalBuses : size_t
    {
      bus, ///< Unique ID of the connecting bus
      SIZE
    };

    /// Signal inputs for a Gensal generator model
    enum class GensalSignalInputs : size_t
    {
      pmech, ///< \f$P_m\f$ Unique ID of the signal providing mechanical power
      efd,   ///< \f$E_{fd}\f$ Unique ID of the signal providing exciter field voltage
      SIZE
    };

    /// Signal outputs for a Gensal generator model
    enum class GensalSignalOutputs : size_t
    {
      speed, ///< \f$\omega\f$ Unique ID of the signal receiving speed deviation
      SIZE
    };

    /// Variables able to be monitored for a Gensal generator model
    enum class GensalMonitorableVariables
    {
      ir,     ///< \f$I_r\f$ Network-frame real terminal current
      ii,     ///< \f$I_i\f$ Network-frame imaginary terminal current
      p,      ///< \f$P\f$ Active power
      q,      ///< \f$Q\f$ Reactive power
      delta,  ///< \f$\delta\f$ Rotor angle
      omega,  ///< \f$\omega\f$ Speed deviation
      speed,  ///< \f$1+\omega\f$ Per-unit machine speed
      Eqp,    ///< \f$E'_q\f$ Q-axis transient voltage
      psidp,  ///< \f$\psi'_d\f$ D-axis transient flux
      psiqpp, ///< \f$\psi''_q\f$ Q-axis subtransient flux
      psidpp, ///< \f$\psi''_d\f$ D-axis subtransient flux
      vd,     ///< \f$V_d\f$ D-axis terminal voltage
      vq,     ///< \f$V_q\f$ Q-axis terminal voltage
      te,     ///< \f$T_e\f$ Electrical torque
      id,     ///< \f$I_d\f$ D-axis current
      iq      ///< \f$I_q\f$ Q-axis current
    };

    /**
     * @brief Contains modeling data for a Gensal generator model.
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     */
    template <typename real_type, typename index_type>
    struct GensalData : public ComponentData<real_type,
                                             index_type,
                                             GensalParameters,
                                             GensalBuses,
                                             GensalSignalInputs,
                                             GensalSignalOutputs,
                                             GensalMonitorableVariables>
    {
      GensalData() = default;

      using Parameters           = GensalParameters;
      using Buses                = GensalBuses;
      using SignalInputs         = GensalSignalInputs;
      using SignalOutputs        = GensalSignalOutputs;
      using MonitorableVariables = GensalMonitorableVariables;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
