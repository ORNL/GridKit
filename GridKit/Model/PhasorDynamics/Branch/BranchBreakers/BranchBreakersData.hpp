/**
 * @file BranchBreakersData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for the breaker-terminated branch model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Parameter keys for the breaker-terminated branch model. Every
    /// parameter is optional and retains its documented default when omitted.
    enum class BranchBreakersParameters : size_t
    {
      R,     ///< \f$R\f$ Branch series resistance [p.u.]
      X,     ///< \f$X\f$ Branch series reactance [p.u.]
      G,     ///< \f$G\f$ Total line shunt conductance, split equally between the two terminals [p.u.]
      B,     ///< \f$B\f$ Total line shunt susceptance, split equally between the two terminals [p.u.]
      Gmag,  ///< \f$G_{\mathrm{mag}}\f$ Magnetizing shunt conductance at bus 1, the tapped side [p.u.]
      Bmag,  ///< \f$B_{\mathrm{mag}}\f$ Magnetizing shunt susceptance at bus 1, the tapped side [p.u.]
      tap,   ///< \f$\tau\f$ Off-nominal tap magnitude on the bus-1 side [p.u.]
      phase, ///< \f$\theta\f$ Off-nominal phase-shift angle [rad]
      Tbrk,  ///< \f$T_{\mathrm{brk}}\f$ Breaker operating time, command to half travel [sec]
    };

    /// Buses for the breaker-terminated branch model.
    enum class BranchBreakersBuses : size_t
    {
      bus1, ///< \f$V_{\mathrm{r}1},V_{\mathrm{i}1}\f$ Required Known bus-1 terminal voltage, the tapped side [p.u.]
      bus2, ///< \f$V_{\mathrm{r}2},V_{\mathrm{i}2}\f$ Required Known bus-2 terminal voltage [p.u.]
    };

    /// Signal inputs for the breaker-terminated branch model.
    enum class BranchBreakersSignalInputs : size_t
    {
      trip1,  ///< \f$s_1\f$ Optional Known bus-1 breaker trip command [-]
      reset1, ///< \f$r_1\f$ Optional Known bus-1 breaker reset command [-]
      trip2,  ///< \f$s_2\f$ Optional Known bus-2 breaker trip command [-]
      reset2, ///< \f$r_2\f$ Optional Known bus-2 breaker reset command [-]
    };

    /// Signal outputs for the breaker-terminated branch model.
    enum class BranchBreakersSignalOutputs : size_t
    {
      ir1, ///< \f$I_{\mathrm{r}1}\f$ Optional Known bus-1 terminal-current real-component output [p.u.]
      ii1, ///< \f$I_{\mathrm{i}1}\f$ Optional Known bus-1 terminal-current imaginary-component output [p.u.]
      ir2, ///< \f$I_{\mathrm{r}2}\f$ Optional Known bus-2 terminal-current real-component output [p.u.]
      ii2, ///< \f$I_{\mathrm{i}2}\f$ Optional Known bus-2 terminal-current imaginary-component output [p.u.]
    };

    /// Variables available through the monitor interface.
    enum class BranchBreakersMonitorableVariables : size_t
    {
      ir1, ///< \f$I_{\mathrm{r}1}\f$ Bus-1 terminal-current real component [p.u.]
      ii1, ///< \f$I_{\mathrm{i}1}\f$ Bus-1 terminal-current imaginary component [p.u.]
      im1, ///< \f$I_{\mathrm{m}1}\f$ Bus-1 terminal-current magnitude [p.u.]
      p1,  ///< \f$P_1\f$ Bus-1 terminal active power [p.u.]
      q1,  ///< \f$Q_1\f$ Bus-1 terminal reactive power [p.u.]
      ir2, ///< \f$I_{\mathrm{r}2}\f$ Bus-2 terminal-current real component [p.u.]
      ii2, ///< \f$I_{\mathrm{i}2}\f$ Bus-2 terminal-current imaginary component [p.u.]
      im2, ///< \f$I_{\mathrm{m}2}\f$ Bus-2 terminal-current magnitude [p.u.]
      p2,  ///< \f$P_2\f$ Bus-2 terminal active power [p.u.]
      q2,  ///< \f$Q_2\f$ Bus-2 terminal reactive power [p.u.]
    };

    /**
     * @brief Model data for breaker-terminated branch parameters, terminal
     *        buses, signal ports, and monitored variables.
     *
     * @tparam real_type Real parameter value type.
     * @tparam index_type Integer index type.
     *
     * @see BranchBreakers
     */
    template <typename real_type, typename index_type>
    using BranchBreakersData =
        ComponentData<real_type,
                      index_type,
                      BranchBreakersParameters,
                      BranchBreakersBuses,
                      BranchBreakersSignalInputs,
                      BranchBreakersSignalOutputs,
                      BranchBreakersMonitorableVariables>;
  } // namespace PhasorDynamics
} // namespace GridKit
