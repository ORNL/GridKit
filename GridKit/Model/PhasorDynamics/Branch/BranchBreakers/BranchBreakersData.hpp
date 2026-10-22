/**
 * @file BranchBreakersData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for the BranchBreakers model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Parameter keys for a BranchBreakers. Every parameter is optional and
    /// retains its documented default when omitted.
    enum class BranchBreakersParameters : size_t
    {
      R,     ///< \f$R\f$ Branch series resistance [p.u.]
      X,     ///< \f$X\f$ Branch series reactance [p.u.]
      G,     ///< \f$G\f$ Total line shunt conductance [p.u.]
      B,     ///< \f$B\f$ Total line shunt susceptance [p.u.]
      Gmag,  ///< \f$G_{\mathrm{mag}}\f$ Magnetizing shunt conductance at bus 1 [p.u.]
      Bmag,  ///< \f$B_{\mathrm{mag}}\f$ Magnetizing shunt susceptance at bus 1 [p.u.]
      tap,   ///< \f$\tau\f$ Off-nominal tap magnitude at bus 1 [p.u.]
      phase, ///< \f$\theta\f$ Phase-shift angle [rad]
      Tbrk,  ///< \f$T_{\mathrm{brk}}\f$ Breaker operating time [sec]
    };

    /// Buses for a BranchBreakers.
    enum class BranchBreakersBuses : size_t
    {
      bus1, ///< \f$V_{\mathrm{r}1},V_{\mathrm{i}1}\f$ Required Known bus-1 terminal voltage, the tapped side [p.u.]
      bus2, ///< \f$V_{\mathrm{r}2},V_{\mathrm{i}2}\f$ Required Known bus-2 terminal voltage [p.u.]
    };

    /// Signal inputs for a BranchBreakers.
    enum class BranchBreakersSignalInputs : size_t
    {
      trip1,  ///< \f$s_1\f$ Optional Known bus-1 breaker trip command [-]
      close1, ///< \f$c_1\f$ Optional Known bus-1 breaker close command [-]
      trip2,  ///< \f$s_2\f$ Optional Known bus-2 breaker trip command [-]
      close2, ///< \f$c_2\f$ Optional Known bus-2 breaker close command [-]
    };

    /// Signal outputs for a BranchBreakers, in computed-output order.
    enum class BranchBreakersSignalOutputs : size_t
    {
      ir1, ///< \f$I_{\mathrm{r}1}\f$ Optional Known bus-1 terminal current, real component [p.u.]
      ii1, ///< \f$I_{\mathrm{i}1}\f$ Optional Known bus-1 terminal current, imaginary component [p.u.]
      ir2, ///< \f$I_{\mathrm{r}2}\f$ Optional Known bus-2 terminal current, real component [p.u.]
      ii2, ///< \f$I_{\mathrm{i}2}\f$ Optional Known bus-2 terminal current, imaginary component [p.u.]
    };

    /// Variables available through the monitor interface.
    enum class BranchBreakersMonitorableVariables : size_t
    {
      ir1, ///< \f$I_{\mathrm{r}1}\f$ Bus-1 terminal current, real component [p.u.]
      ii1, ///< \f$I_{\mathrm{i}1}\f$ Bus-1 terminal current, imaginary component [p.u.]
      im1, ///< \f$I_{\mathrm{m}1}\f$ Bus-1 terminal current magnitude [p.u.]
      p1,  ///< \f$P_1\f$ Active power entering bus 1 [p.u.]
      q1,  ///< \f$Q_1\f$ Reactive power entering bus 1 [p.u.]
      ir2, ///< \f$I_{\mathrm{r}2}\f$ Bus-2 terminal current, real component [p.u.]
      ii2, ///< \f$I_{\mathrm{i}2}\f$ Bus-2 terminal current, imaginary component [p.u.]
      im2, ///< \f$I_{\mathrm{m}2}\f$ Bus-2 terminal current magnitude [p.u.]
      p2,  ///< \f$P_2\f$ Active power entering bus 2 [p.u.]
      q2,  ///< \f$Q_2\f$ Reactive power entering bus 2 [p.u.]
    };

    /**
     * @brief Model data for a BranchBreakers: parameter values, terminal
     *        buses, breaker commands, terminal-current outputs, and
     *        monitored-variable selections.
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
