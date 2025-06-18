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
      /// Load resistance
      RealT R{0.0};

      /// Load reactance
      RealT X{0.0};

      //TODO: add the monitorable variables to this

      /// Unique ID of bus to which the load is connnected
      IdxT bus_id{0};
    };
  } // namespace PhasorDynamics
} // namespace GridKit
