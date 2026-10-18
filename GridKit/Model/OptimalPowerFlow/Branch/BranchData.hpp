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
      R,          ///< \f$R\f$ Series resistance [p.u.]
      X,          ///< \f$X\f$ Series reactance [p.u.]
      G,          ///< \f$G\f$ Total line shunt conductance, split between the ports [p.u.]
      B,          ///< \f$B\f$ Total line shunt susceptance, split between the ports [p.u.]
      Gmag,       ///< \f$G_\mathrm{mag}\f$ Magnetizing shunt conductance at bus 1 [p.u.]
      Bmag,       ///< \f$B_\mathrm{mag}\f$ Magnetizing shunt susceptance at bus 1 [p.u.]
      tap,        ///< \f$\tau\f$ Off-nominal tap magnitude on the bus-1 side [p.u.]
      phase,      ///< \f$\theta\f$ Phase-shift angle [rad]
      Smax,       ///< \f$S^{\max}\f$ Apparent power limit at each port [p.u.]
      in_service, ///< \f$u\f$ Service status
    };

    /// Ports of a branch
    enum class BranchPorts : size_t
    {
      bus1, ///< Bus-1 port, the tapped side
      bus2, ///< Bus-2 port
    };

    /**
     * @brief Contains modeling data for a Branch
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer data type
     */
    template <typename real_type, typename index_type>
    using BranchData = ComponentData<real_type, index_type, BranchParameters, BranchPorts>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
