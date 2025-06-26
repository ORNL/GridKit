/**
 * @file BusControlData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Modeling data for buses (nodes)
 *
 */
#pragma once

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Contains modeling data for a Bus
     *
     * @tparam RealT Real parameter data type
     * @tparam IdxT  Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     *
     * @todo Decide on naming scheme for model parameters.
     */
    template <typename RealT, typename IdxT>
    struct BusControlData
    {
      enum BusType
      {
        INVALID = 0,
        SIGNAL      ///< Control   , Signal
      };

      IdxT    bus_id{0}; ///< Unique ID of bus 1
      BusType bus_type{INVALID};
    };
  } // namespace PhasorDynamics
} // namespace GridKit
