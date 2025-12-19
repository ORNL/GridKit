/**
 * @file BranchData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Modeling data for branches (transmission lines)
 *
 */
#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Initial parameters for a branch
    enum class BranchParameters
    {
      R, ///< Line series resistance
      X, ///< Line series reactance
      G, ///< Line shunt conductance
      B, ///< Line shunt charging
    };

    /// Ports for a branch
    enum class BranchPorts
    {
      bus1, ///< Unique ID of bus 1
      bus2, ///< Unique ID of bus 2
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
      q2,
      SIZE
    };

    /**
     * @brief Convert enum value to string label
     */
    inline const std::string& enumLabel(BranchMonitorableVariables v)
    {
      static const std::string labels[] = {
          "ir1",
          "ii1",
          "im1",
          "p1",
          "q1",
          "ir2",
          "ii2",
          "im2",
          "p2",
          "q2"};
      return labels[Utilities::enumId(v)];
    }

    /**
     * @brief Contains modeling data for a Branch
     *
     * @tparam RealT Real parameter data type
     * @tparam IdxT  Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     */
    template <typename RealT, typename IdxT>
    struct BranchData : public ComponentData<RealT,
                                             IdxT,
                                             BranchParameters,
                                             BranchPorts,
                                             BranchMonitorableVariables>
    {
      BranchData() = default;

      using Parameters           = BranchParameters;
      using Ports                = BranchPorts;
      using MonitorableVariables = BranchMonitorableVariables;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
