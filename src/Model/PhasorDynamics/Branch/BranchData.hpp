/**
 * @file BranchData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Modeling data for branches (transmission lines)
 *
 */
#pragma once

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Contains modeling data for a Branch
     *
     * @tparam RealT Real parameter data type
     * @tparam IdxT  Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     *
     * @todo Decide on naming scheme for model parameters.
     */
    template <typename RealT, typename IdxT>
    struct BranchData
    {
      RealT R{0.0}; ///< line series resistance
      RealT X{0.0}; ///< line series reactance
      RealT G{0.0}; ///< line shunt conductance
      RealT B{0.0}; ///< line shunt charging

      IdxT bus1_id{0}; ///< Unique ID of bus 1
      IdxT bus2_id{0}; ///< Unique ID of bus 2
    };
  } // namespace PhasorDynamics
} // namespace GridKit
