/**
 * @file ComponentData.hpp
 * @brief Data of optimal power flow components.
 */

#pragma once

#include <array>
#include <cstddef>
#include <map>
#include <string>

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

      /// An out-of-service component keeps its entries and contributes nothing
      bool in_service{true};

      /// Mapping of parameters to parameter values
      std::map<Parameters, RealT> parameters;

      /// Mapping of ports to the numbers of the components they connect to
      std::map<Ports, IdxT> ports;
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
     * @brief Number at each port in `Ports` order
     *
     * A port without one gets `INVALID_INDEX`.
     */
    template <typename DataT>
    auto portNumbers(const DataT& data)
    {
      using IdxT = typename DataT::IdxT;

      std::array<IdxT, Utilities::enum_size<typename DataT::Ports>()> numbers;
      numbers.fill(INVALID_INDEX<IdxT>);
      for (const auto& [port, number] : data.ports)
      {
        numbers[static_cast<size_t>(port)] = number;
      }
      return numbers;
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
