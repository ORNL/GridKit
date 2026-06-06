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
    enum class BranchParameters
    {
      R,     ///< Line series resistance
      X,     ///< Line series reactance
      G,     ///< Total shunt conductance
      B,     ///< Total shunt susceptance
      tap,   ///< Off-nominal tap magnitude on bus1 side
      phase, ///< Phase shift angle in radians
    };

    /// Terminals for a branch
    enum class BranchTerminals : size_t
    {
      bus1, ///< Unique ID of bus 1
      bus2, ///< Unique ID of bus 2
      SIZE
    };

    /// Input ports supported for a branch
    enum class BranchInputPorts : size_t
    {
      SIZE
    };

    /// Output ports supported for a branch
    enum class BranchOutputPorts : size_t
    {
      SIZE
    };

    /// Variables able to be monitored for a branch
    enum class BranchMonitorableVariables
    {
      ir1,
      ii1,
      im1,
      p1,
      q1,
      ir2,
      ii2,
      im2,
      p2,
      q2
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
    struct BranchData : public ComponentData<real_type,
                                             index_type,
                                             BranchParameters,
                                             BranchTerminals,
                                             BranchInputPorts,
                                             BranchOutputPorts,
                                             BranchMonitorableVariables>
    {
      BranchData() = default;

      using Parameters           = BranchParameters;
      using Terminals            = BranchTerminals;
      using InputPorts           = BranchInputPorts;
      using OutputPorts          = BranchOutputPorts;
      using MonitorableVariables = BranchMonitorableVariables;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
