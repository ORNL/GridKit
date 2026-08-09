/**
 * @file ForcedOscillationStudy.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Study-time forced-oscillation injection utilities.
 */

#pragma once

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <magic_enum/magic_enum.hpp>
#include <nlohmann/json.hpp>

#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief One forced-oscillation injection requested by a study file.
     *
     * The target uses the form `Class:id.port`, where `port` names an entry
     * in the target component data type's `SignalInputs` enum. The only
     * supported mode is `add`.
     */
    struct ForcedOscillationInjection
    {
      std::string              id;                                  ///< Unique forced-oscillation component ID
      std::string              target;                              ///< Consumer target in `Class:id.port` form
      std::string              mode{"add"};                         ///< Injection mode; currently only `add`
      nlohmann::json           params{nlohmann::json::object()};    ///< Source parameters
      std::vector<std::string> mon{"output", "envelope", "active"}; ///< Monitors
    };

    /**
     * @brief Parse one forced-oscillation study entry.
     */
    inline void from_json(const nlohmann::json& j, ForcedOscillationInjection& injection)
    {
      injection = ForcedOscillationInjection{};
      j.at("id").get_to(injection.id);
      j.at("target").get_to(injection.target);
      injection.mode = j.value("mode", std::string{"add"});

      if (j.contains("params"))
      {
        injection.params = j.at("params");
        if (!injection.params.is_object())
        {
          throw std::runtime_error("Forced oscillation params must be an object");
        }
      }

      if (j.contains("mon"))
      {
        j.at("mon").get_to(injection.mon);
      }
    }

    namespace detail
    {
      [[noreturn]] inline void forcedOscillationStudyError(const std::string& message)
      {
        ::GridKit::Utilities::Logger::error() << message << '\n';
        throw std::runtime_error(message);
      }

      struct ForcedOscillationTarget
      {
        std::string component_class;
        std::string component_id;
        std::string signal_input;
      };

      inline ForcedOscillationTarget parseForcedOscillationTarget(const std::string& target)
      {
        const auto colon = target.find(':');
        const auto dot   = target.rfind('.');
        if (colon == std::string::npos || dot == std::string::npos
            || colon == 0 || dot <= colon + 1 || dot + 1 == target.size())
        {
          forcedOscillationStudyError(
              "Malformed forced oscillation target \"" + target
              + "\"; expected Class:id.port");
        }

        return {
            target.substr(0, colon),
            target.substr(colon + 1, dot - colon - 1),
            target.substr(dot + 1)};
      }

      template <typename DataT>
      DataT* findForcedOscillationTarget(
          std::vector<DataT>& devices,
          const std::string&  id)
      {
        auto device = std::find_if(
            devices.begin(),
            devices.end(),
            [&id](const DataT& candidate)
            { return candidate.disambiguation_string == id; });
        if (device == devices.end())
        {
          return nullptr;
        }

        const auto duplicate = std::find_if(
            std::next(device),
            devices.end(),
            [&id](const DataT& candidate)
            { return candidate.disambiguation_string == id; });
        if (duplicate != devices.end())
        {
          forcedOscillationStudyError(
              "Ambiguous forced oscillation target id \"" + id + "\"");
        }
        return &(*device);
      }

      template <typename ModelDataT>
      std::set<typename ModelDataT::IdxT> collectSignalIds(const ModelDataT& model)
      {
        std::set<typename ModelDataT::IdxT> signal_ids;
        for (const auto& signal : model.signal)
        {
          if (!signal_ids.insert(signal.signal_id).second)
          {
            forcedOscillationStudyError(
                "Duplicate signal ID " + std::to_string(signal.signal_id));
          }
        }
        return signal_ids;
      }

      template <typename ModelDataT>
      typename ModelDataT::SignalDataT& appendSignal(
          ModelDataT&                          model,
          std::set<typename ModelDataT::IdxT>& signal_ids,
          const std::string&                   name)
      {
        using IdxT = typename ModelDataT::IdxT;

        IdxT signal_id{0};
        if (!signal_ids.empty())
        {
          const IdxT maximum_id = *signal_ids.rbegin();
          if (maximum_id == std::numeric_limits<IdxT>::max())
          {
            forcedOscillationStudyError(
                "Unable to allocate a signal ID for forced oscillation \""
                + name + "\"");
          }
          signal_id = maximum_id + IdxT{1};
        }

        auto& signal     = model.signal.emplace_back();
        signal.signal_id = signal_id;
        signal.name      = name;
        signal_ids.insert(signal_id);
        return signal;
      }

      template <typename ModelDataT>
      void appendForcedOscillationSource(
          ModelDataT&                       model,
          const ForcedOscillationInjection& injection,
          typename ModelDataT::IdxT         output_signal)
      {
        using DataT         = typename ModelDataT::ForcedOscillationDataT;
        using RealT         = typename ModelDataT::RealT;
        using Parameters    = typename DataT::Parameters;
        using SignalOutputs = typename DataT::SignalOutputs;
        using Monitors      = typename DataT::MonitorableVariables;

        DataT source;
        source.device_class                          = "ForcedOscillation";
        source.disambiguation_string                 = injection.id;
        source.signal_outputs[SignalOutputs::output] = output_signal;

        for (const auto& [name, value] : injection.params.items())
        {
          const auto parameter = magic_enum::enum_cast<Parameters>(name);
          if (!parameter.has_value())
          {
            forcedOscillationStudyError(
                "Unknown ForcedOscillation parameter \"" + name
                + "\" in injection \"" + injection.id + "\"");
          }
          if (!value.is_number())
          {
            forcedOscillationStudyError(
                "ForcedOscillation parameter \"" + name
                + "\" in injection \"" + injection.id
                + "\" must be numeric");
          }
          source.parameters[*parameter] = value.template get<RealT>();
        }

        for (const auto& name : injection.mon)
        {
          if (name == "output")
          {
            source.monitored_variables.insert(Monitors::output);
          }
          else if (name == "envelope")
          {
            source.monitored_variables.insert(Monitors::envelope);
          }
          else if (name == "active")
          {
            source.monitored_variables.insert(Monitors::active);
          }
          else
          {
            forcedOscillationStudyError(
                "Unknown ForcedOscillation monitor \"" + name
                + "\" in injection \"" + injection.id + "\"");
          }
        }

        model.forced_oscillation.push_back(std::move(source));
      }

      template <typename ModelDataT, typename DataT>
      void injectAtSignalInput(
          ModelDataT&                          model,
          std::vector<DataT>&                  devices,
          const ForcedOscillationTarget&       target,
          const ForcedOscillationInjection&    injection,
          std::set<typename ModelDataT::IdxT>& signal_ids)
      {
        using RealT        = typename ModelDataT::RealT;
        using SignalInputs = typename DataT::SignalInputs;

        const auto signal_input = magic_enum::enum_cast<SignalInputs>(target.signal_input);
        if (!signal_input.has_value() || *signal_input == SignalInputs::SIZE)
        {
          forcedOscillationStudyError(
              "Unknown signal input \"" + target.signal_input
              + "\" on target \"" + injection.target + "\"");
        }

        auto* device = findForcedOscillationTarget(devices, target.component_id);
        if (device == nullptr)
        {
          forcedOscillationStudyError(
              "Unable to find " + target.component_class + " device \""
              + target.component_id + "\"");
        }
        if (!device->signal_inputs.contains(*signal_input))
        {
          forcedOscillationStudyError(
              "add mode requires an existing signal on target \""
              + injection.target + "\"");
        }

        const auto original_signal = device->signal_inputs.at(*signal_input);
        if (!signal_ids.contains(original_signal))
        {
          forcedOscillationStudyError(
              "Target \"" + injection.target + "\" references undeclared signal "
              + std::to_string(original_signal));
        }

        auto& forcing = appendSignal(
            model, signal_ids, injection.id + "_forcing");
        const auto forcing_signal = forcing.signal_id;

        auto& junction = appendSignal(
            model, signal_ids, injection.id + "_junction");
        const auto junction_signal         = junction.signal_id;
        auto&      junction_data           = junction.junction.emplace();
        junction_data.bias                 = RealT{0};
        junction_data.initialization_input = original_signal;
        junction_data.inputs.push_back({original_signal, RealT{1}});
        junction_data.inputs.push_back({forcing_signal, RealT{1}});

        appendForcedOscillationSource(model, injection, forcing_signal);
        device->signal_inputs[*signal_input] = junction_signal;
      }

      template <typename ModelDataT>
      void dispatchForcedOscillationTarget(
          ModelDataT&                          model,
          const ForcedOscillationTarget&       target,
          const ForcedOscillationInjection&    injection,
          std::set<typename ModelDataT::IdxT>& signal_ids)
      {
        bool matched  = false;
        auto dispatch = [&]<typename DataT>(
                            std::string_view    component_class,
                            std::vector<DataT>& devices)
        {
          if (!matched && target.component_class == component_class)
          {
            injectAtSignalInput(model, devices, target, injection, signal_ids);
            matched = true;
          }
        };

        dispatch("BusToSignalAdapter", model.adapter);
        dispatch("BusFault", model.bus_fault);
        dispatch("Regca", model.regca);
        dispatch("Reecb", model.reecb);
        dispatch("Repca", model.repca);
        dispatch("Genrou", model.genrou);
        dispatch("Gensal", model.gensal);
        dispatch("GenClassical", model.genclassical);
        dispatch("Tgov1", model.gov);
        dispatch("Esdc1a", model.esdc1a);
        dispatch("GastPti", model.gastpti);
        dispatch("Hygov", model.hygov);
        dispatch("Ieeet1", model.exciter);
        dispatch("SexsPti", model.sexspti);
        dispatch("Ieeest", model.stabilizer);

        if (!matched)
        {
          forcedOscillationStudyError(
              "Unknown forced oscillation target class \""
              + target.component_class + "\"");
        }
      }

    } // namespace detail

    /**
     * @brief Apply additive forced oscillations to a copy of model data.
     *
     * Each injection adds one source signal and one algebraic junction signal.
     * The junction combines the target's original input with the source, and
     * the typed target `signal_inputs` entry is redirected to that junction.
     */
    template <typename RealT, typename IdxT>
    SystemModelData<RealT, IdxT> applyForcedOscillations(
        const SystemModelData<RealT, IdxT>&            base,
        const std::vector<ForcedOscillationInjection>& injections)
    {
      auto model      = base;
      auto signal_ids = detail::collectSignalIds(model);

      std::set<std::string> injection_ids;
      std::set<std::string> injection_targets;
      for (const auto& source : model.forced_oscillation)
      {
        if (!injection_ids.insert(source.disambiguation_string).second)
        {
          detail::forcedOscillationStudyError(
              "Duplicate ForcedOscillation component id \""
              + source.disambiguation_string + "\"");
        }
      }

      for (const auto& injection : injections)
      {
        if (injection.id.empty())
        {
          detail::forcedOscillationStudyError(
              "Forced oscillation injection id cannot be empty");
        }
        if (!injection_ids.insert(injection.id).second)
        {
          detail::forcedOscillationStudyError(
              "Duplicate forced oscillation injection id \""
              + injection.id + "\"");
        }
        if (injection.mode != "add")
        {
          detail::forcedOscillationStudyError(
              "Unknown forced oscillation mode \"" + injection.mode
              + "\"; only \"add\" is supported");
        }
        if (!injection_targets.insert(injection.target).second)
        {
          detail::forcedOscillationStudyError(
              "Duplicate forced oscillation target \""
              + injection.target + "\"");
        }

        const auto target = detail::parseForcedOscillationTarget(injection.target);
        detail::dispatchForcedOscillationTarget(
            model, target, injection, signal_ids);
      }

      return model;
    }
  } // namespace PhasorDynamics
} // namespace GridKit
