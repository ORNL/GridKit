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

    /// Data of a device: its `id`, parameter values, and the bus number or signal id at each port
    template <typename real_type,
              typename index_type,
              typename parameters_enum,
              typename ports_enum,
              typename signals_enum = Utilities::EmptyEnum>
      requires Utilities::SizedEnum<parameters_enum>
               && Utilities::SizedEnum<ports_enum>
               && Utilities::SizedEnum<signals_enum>
    struct ComponentData
    {
      using RealT      = real_type;
      using IdxT       = index_type;
      using Parameters = parameters_enum;
      using Ports      = ports_enum;
      using Signals    = signals_enum;

      std::string id;

      ParameterMap<Parameters, RealT, IdxT> parameters;

      std::map<Ports, IdxT> ports;

      std::map<Signals, IdxT> signals;
    };
  } // namespace OptimalPowerFlow
} // namespace GridKit
