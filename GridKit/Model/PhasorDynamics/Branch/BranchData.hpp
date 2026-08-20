/**
 * @file BranchData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Modeling data for branches
 *
 */
#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Initial parameters for a branch
    enum class BranchParameters : size_t
    {
      R,     ///< \f$R\f$ Branch series resistance [p.u.]
      X,     ///< \f$X\f$ Branch series reactance [p.u.]
      G,     ///< \f$G\f$ Total line shunt conductance, split equally between the two terminals [p.u.]
      B,     ///< \f$B\f$ Total line shunt susceptance, split equally between the two terminals [p.u.]
      Gmag,  ///< \f$G_{\mathrm{mag}}\f$ Magnetizing shunt conductance at bus 1, the tapped side [p.u.]
      Bmag,  ///< \f$B_{\mathrm{mag}}\f$ Magnetizing shunt susceptance at bus 1, the tapped side [p.u.]
      tap,   ///< \f$\tau\f$ Off-nominal tap magnitude on the bus-1 side [p.u.]
      phase, ///< \f$\theta\f$ Off-nominal phase-shift angle [rad]
    };

    /// Buses for a branch
    enum class BranchBuses : size_t
    {
      bus1, ///< \f$V_{\mathrm{r}1},V_{\mathrm{i}1}\f$ Required Known bus-1 terminal voltage, the tapped side [p.u.]
      bus2, ///< \f$V_{\mathrm{r}2},V_{\mathrm{i}2}\f$ Required Known bus-2 terminal voltage [p.u.]
    };

    /// Signal inputs supported for a branch
    enum class BranchSignalInputs : size_t
    {
    };

    /// Signal outputs supported for a branch
    enum class BranchSignalOutputs : size_t
    {
    };

    /// Variables able to be monitored for a branch
    enum class BranchMonitorableVariables : size_t
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
     * @brief Contains modeling data for a Branch
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     */
    template <typename real_type, typename index_type>
    using BranchData =
        ComponentData<real_type,
                      index_type,
                      BranchParameters,
                      BranchBuses,
                      BranchSignalInputs,
                      BranchSignalOutputs,
                      BranchMonitorableVariables>;
  } // namespace PhasorDynamics
} // namespace GridKit
