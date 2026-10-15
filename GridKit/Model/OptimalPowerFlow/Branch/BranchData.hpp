/**
 * @file BranchData.hpp
 * @brief Modeling data for optimal power flow branches.
 */

#pragma once

#include <cstddef>

#include <GridKit/Model/OptimalPowerFlow/ComponentData.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /// Parameters for a branch
    enum class BranchParameters : size_t
    {
      R,     ///< \f$R\f$ Series resistance [p.u.]
      X,     ///< \f$X\f$ Series reactance [p.u.]
      G,     ///< \f$G\f$ Total line shunt conductance, split between the terminals [p.u.]
      B,     ///< \f$B\f$ Total line shunt susceptance, split between the terminals [p.u.]
      Gmag,  ///< \f$G_\mathrm{mag}\f$ Magnetizing shunt conductance at bus 1 [p.u.]
      Bmag,  ///< \f$B_\mathrm{mag}\f$ Magnetizing shunt susceptance at bus 1 [p.u.]
      tap,   ///< \f$\tau\f$ Off-nominal tap magnitude on the bus-1 side [p.u.]
      phase, ///< \f$\theta\f$ Phase-shift angle [rad]
      Smax,  ///< \f$S^{\max}\f$ Apparent power limit at each terminal [p.u.]
    };

    /// Buses for a branch
    enum class BranchBuses : size_t
    {
      bus1, ///< Bus-1 terminal, the tapped side
      bus2, ///< Bus-2 terminal
    };

    /**
     * @brief Contains modeling data for a Branch
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer data type
     */
    template <typename real_type, typename index_type>
    using BranchData = ComponentData<real_type, index_type, BranchParameters, BranchBuses>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
