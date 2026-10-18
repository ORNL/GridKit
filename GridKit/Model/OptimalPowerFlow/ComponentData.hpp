/**
 * @file ComponentData.hpp
 * @brief Data of optimal power flow components.
 */

#pragma once

#include <map>
#include <string>
#include <variant>

#include <GridKit/Utilities/Enum.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /// Value of each parameter, typed as in the input file
    template <typename Parameters, typename RealT, typename IdxT>
    using ParameterMap = std::map<Parameters, std::variant<std::string, bool, RealT, IdxT>>;

    /**
     * @brief Unified interface for `Component` data containers
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer data type
     * @tparam parameters_enum Enumeration over parameters
     * @tparam ports_enum Enumeration over ports
     */
    template <typename real_type, typename index_type, typename parameters_enum, typename ports_enum>
      requires Utilities::SizedEnum<parameters_enum> && Utilities::SizedEnum<ports_enum>
    struct ComponentData
    {
      using RealT      = real_type;
      using IdxT       = index_type;
      using Parameters = parameters_enum;
      using Ports      = ports_enum;

      /// Component `id`
      std::string id;

      /// Mapping of parameters to parameter values
      ParameterMap<Parameters, RealT, IdxT> parameters;

      /// Mapping of ports to the numbers of the buses they connect to
      std::map<Ports, IdxT> ports;
    };
  } // namespace OptimalPowerFlow
} // namespace GridKit
