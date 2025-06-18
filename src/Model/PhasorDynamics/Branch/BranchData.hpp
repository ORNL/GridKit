/**
 * @file BranchData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Modeling data for branches (transmission lines)
 *
 */
#pragma once

#include <bitset>
#include <type_traits>

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
      /// Line series resistance
      RealT R{0.0};

      /// Line series reactance
      RealT X{0.0};

      /// Line shunt conductance
      RealT G{0.0};

      /// Line shunt charging
      RealT B{0.0};

      /// Unique ID of bus 1
      IdxT bus1_id{0};

      /// Unique ID of bus 2
      IdxT bus2_id{0};

      /// Indices of the variables able to be monitored in the bitset
      enum class MonitorableVariables : size_t
      {
        Vm,
        Va,
        Maximum,
      };

      /// Set of variables being monitored
      std::bitset<static_cast<std::underlying_type_t<MonitorableVariables>>(MonitorableVariables::Maximum) - 1> monitored_variables;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
