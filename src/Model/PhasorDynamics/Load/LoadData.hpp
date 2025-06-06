/**
 * @file LoadData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Modeling data for loads
 *
 */
#pragma once

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Contains modeling data for a Load
     *
     * @tparam RealT Real parameter data type
     * @tparam IdxT  Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     *
     * @todo Decide on naming scheme for model parameters.
     */
    template <typename RealT, typename IdxT>
    struct LoadData
    {
      RealT R{0.0}; ///< load resistance
      RealT X{0.0}; ///< load reactance

      IdxT bus_id{0}; ///< Unique ID of bus to which the load is connnected.
    };
  } // namespace PhasorDynamics
} // namespace GridKit
