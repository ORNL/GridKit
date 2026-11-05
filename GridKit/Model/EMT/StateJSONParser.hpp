/**
 * @file StateJSONParser.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Initial values of an EMT system from a state file.
 */

#pragma once

#include <algorithm>
#include <format>
#include <stdexcept>
#include <vector>

#include <nlohmann/json.hpp>

#include <GridKit/Model/EMT/ComponentLibrary.hpp>
#include <GridKit/Model/EMT/SystemModel.hpp>
#include <GridKit/ScalarTraits.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit
{
  namespace EMT
  {
    using json = nlohmann::json;

    /// Internal variables and submodels of `model` by enumerator name
    template <typename RealT, typename ModelT>
    void applyEntries(const json& entries, ModelT& model)
    {
      for (const auto& [name, value] : entries.items())
      {
        if constexpr (requires { typename ModelT::SubmodelsT; })
        {
          if (value.is_object())
          {
            applyEntries<RealT>(value, model.submodel(Utilities::enumKey<typename ModelT::SubmodelsT>(name, "submodel", model.id())));
            continue;
          }
        }

        const auto y      = model.y(Utilities::enumKey<typename ModelT::InternalVariablesT>(name, "variable", model.id()));
        const auto values = value.template get<std::vector<RealT>>();
        if (values.size() != y.size())
        {
          throw std::invalid_argument(std::format("{}: '{}' needs {} values", model.id(), name, y.size()));
        }
        std::ranges::copy(values, y.begin());
      }
    }

    /// Buses by number and devices by id; omitted values stay zero
    template <typename ScalarT, typename IdxT>
    void applyState(const json& state, SystemModel<ScalarT, IdxT>& system)
    {
      using RealT = typename ScalarTraits<ScalarT>::RealT;

      const json buses   = state.value("buses", json::object());
      const json devices = state.value("devices", json::object());
      for (const auto& bus : buses.items())
      {
        ComponentLibrary<ScalarT, IdxT>::visit(system.bus(bus.key()), [&bus](auto& model)
                                               { applyEntries<RealT>(bus.value(), model); });
      }
      for (const auto& device : devices.items())
      {
        ComponentLibrary<ScalarT, IdxT>::visit(system.device(device.key()), [&device](auto& model)
                                               { applyEntries<RealT>(device.value(), model); });
      }
    }
  } // namespace EMT
} // namespace GridKit
