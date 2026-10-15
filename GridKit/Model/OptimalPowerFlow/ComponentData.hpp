/**
 * @file ComponentData.hpp
 * @brief Data of optimal power flow components.
 */

#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include <GridKit/Constants.hpp>
#include <GridKit/Utilities/Enum.hpp>

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
      requires Utilities::SizedEnum<parameters_enum> && Utilities::SizedEnum<buses_enum>
    struct ComponentData
    {
      using RealT      = real_type;
      using IdxT       = index_type;
      using Parameters = parameters_enum;
      using Buses      = buses_enum;

      /// Device `id`
      std::string id;

      /// An out-of-service device keeps its entries and carries no power
      bool in_service{true};

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

      std::vector<IdxT> numbers(Utilities::enum_size<typename DataT::Buses>(), INVALID_INDEX<IdxT>);
      for (const auto& [terminal, number] : data.buses)
      {
        numbers[static_cast<size_t>(terminal)] = number;
      }
      return numbers;
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
