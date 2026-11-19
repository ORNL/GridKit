/**
 * @file StateJSONParser.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Initial values of an EMT system from a state file.
 */

#pragma once

#include <algorithm>
#include <format>
#include <span>
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

    /// Internal variables, held inputs and submodels of `model` by enumerator name
    template <typename ScalarT, typename IdxT, typename ModelT>
    void applyEntries(const json& entries, SystemModel<ScalarT, IdxT>& system, ModelT& model)
    {
      using RealT = typename ScalarTraits<ScalarT>::RealT;

      static_assert(Utilities::enum_names_distinct<typename ModelT::InputsT, typename ModelT::InternalVariablesT>(),
                    "an input cannot share its name with an internal variable");

      for (const auto& [name, value] : entries.items())
      {
        if constexpr (requires { typename ModelT::SubmodelsT; })
        {
          if (value.is_object())
          {
            applyEntries(value, system, model.submodel(Utilities::enumKey<typename ModelT::SubmodelsT>(name, "submodel", model.id())));
            continue;
          }
        }

        const auto values = value.template get<std::vector<RealT>>();
        const auto input  = Utilities::enum_cast<typename ModelT::InputsT>(name);
        if (input.has_value())
        {
          system.set(model, input.value(), values);
          continue;
        }

        const auto y = model.y(Utilities::enumKey<typename ModelT::InternalVariablesT>(name, "variable", model.id()));
        if (values.size() != y.size())
        {
          throw std::invalid_argument(std::format("{}: '{}' needs {} values", model.id(), name, y.size()));
        }
        std::ranges::copy(values, y.begin());
      }
    }

    /// Buses by number and devices by id; omitted values are unchanged
    template <typename ScalarT, typename IdxT>
    void applyState(const json& state, SystemModel<ScalarT, IdxT>& system)
    {
      const json buses   = state.value("buses", json::object());
      const json devices = state.value("devices", json::object());
      for (const auto& bus : buses.items())
      {
        ComponentLibrary<ScalarT, IdxT>::visit(system.bus(bus.key()), [&](auto& model)
                                               { applyEntries(bus.value(), system, model); });
      }
      for (const auto& device : devices.items())
      {
        ComponentLibrary<ScalarT, IdxT>::visit(system.device(device.key()), [&](auto& model)
                                               { applyEntries(device.value(), system, model); });
      }
    }
  } // namespace EMT
} // namespace GridKit
