/**
 * @file BusElectricData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Modeling data for buses (electrical)
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
    struct BusElectricData
    {
      enum BusType
      {
        INVALID = 0,
        DEFAULT,    ///< Electrical, Network
        SLACK,      ///< Electrical, Infinite
      };

      RealT Vr0{0.0}; ///< Initial value for real bus voltage
      RealT Vi0{0.0}; ///< Initial value for imaginary bus voltage

      IdxT    bus_id{0}; ///< Unique ID of bus 1
      BusType bus_type{INVALID};
    };
  } // namespace PhasorDynamics
} // namespace GridKit
