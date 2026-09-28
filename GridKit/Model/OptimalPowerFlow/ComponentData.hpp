/**
 * @file ComponentData.hpp
 * @brief Data of optimal power flow components.
 */

#pragma once

#include <map>
#include <string>
#include <type_traits>
#include <vector>

#include <magic_enum/magic_enum.hpp>

#include <GridKit/Constants.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /**
     * @brief Unified interface for `Component` data containers
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer data type
     * @tparam parameters_enum Enumeration over parameters
     * @tparam buses_enum Enumeration over terminal buses
     */
    template <typename real_type, typename index_type, typename parameters_enum, typename buses_enum>
      requires std::is_enum_v<parameters_enum> && std::is_enum_v<buses_enum>
    struct ComponentData
    {
      using RealT      = real_type;
      using IdxT       = index_type;
      using Parameters = parameters_enum;
      using Buses      = buses_enum;

      /// Case device `id`
      std::string id;

      /// Mapping of parameters to parameter values
      std::map<Parameters, RealT> parameters;

      /// Mapping of terminals to bus numbers
      std::map<Buses, IdxT> buses;
    };

    /**
     * @brief Value of parameter `key`, or `fallback` if it is not set
     */
    template <typename DataT>
    typename DataT::RealT parameter(const DataT&               data,
                                    typename DataT::Parameters key,
                                    typename DataT::RealT      fallback)
    {
      const auto entry = data.parameters.find(key);
      if (entry == data.parameters.end())
      {
        return fallback;
      }
      return entry->second;
    }

    /**
     * @brief Bus number of each terminal in `Buses` order
     *
     * A terminal without a bus gets `INVALID_INDEX`.
     */
    template <typename DataT>
    std::vector<typename DataT::IdxT> busNumbers(const DataT& data)
    {
      using IdxT = typename DataT::IdxT;

      std::vector<IdxT> numbers;
      for (const auto terminal : magic_enum::enum_values<typename DataT::Buses>())
      {
        IdxT       number = INVALID_INDEX<IdxT>;
        const auto entry  = data.buses.find(terminal);
        if (entry != data.buses.end())
        {
          number = entry->second;
        }
        numbers.push_back(number);
      }
      return numbers;
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
