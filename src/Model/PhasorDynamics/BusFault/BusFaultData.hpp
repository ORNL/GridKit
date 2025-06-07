/**
 * @file BusFaultData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Modeling data for short-to-ground fault
 *
 */
#pragma once

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Contains modeling data for a short-to-ground fault
     *
     * @tparam RealT Real parameter data type
     * @tparam IdxT  Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     *
     * @todo Decide on naming scheme for model parameters.
     */
    template <typename RealT, typename IdxT>
    struct BusFaultData
    {
      RealT R{0.0}; ///< short to ground resistance
      RealT X{0.0}; ///< short to ground reactance
      int status{0}; ///< if the fault happened

      IdxT bus_id{0}; ///< Unique ID of bus where fault occurs.
    };
  } // namespace PhasorDynamics
} // namespace GridKit
