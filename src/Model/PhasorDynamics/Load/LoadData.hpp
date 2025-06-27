/**
 * @file LoadData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Modeling data for loads
 *
 */
#pragma once

#include <optional>
#include <string>

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
      RealT R{0.0}; ///< Load resistance
      RealT X{0.0}; ///< Load reactance

      IdxT bus_id{0}; ///< Unique ID of bus to which the load is connnected

      std::optional<RealT> freq_base; ///< Override for the system-wide base frequency
      std::optional<RealT> va_base;   ///< Override for the system-wide power base

      std::string disambiguation_string; ///< Disambiguation string for this device

      // TODO: add the monitorable variables to this
    };
  } // namespace PhasorDynamics
} // namespace GridKit
