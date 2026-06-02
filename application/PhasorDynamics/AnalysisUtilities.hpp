#pragma once

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <magic_enum/magic_enum.hpp>
#include <nlohmann/json.hpp>

#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace fs = ::std::filesystem;

    using json = ::nlohmann::json;
    using Log  = GridKit::Utilities::Logger;

    namespace detail
    {
      [[noreturn]] inline void studyError(const std::string& message)
      {
        Log::error() << message << std::endl;
        throw std::runtime_error(message);
      }
    } // namespace detail

    enum class SignalPortRole
    {
      Consumer,
      Producer,
      NonSignal
    };

    inline std::optional<SignalPortRole> portRole(BusToSignalAdapterPorts port)
    {
      switch (port)
      {
      case BusToSignalAdapterPorts::bus:
        return SignalPortRole::NonSignal;
      case BusToSignalAdapterPorts::vr:
      case BusToSignalAdapterPorts::vi:
        return SignalPortRole::Producer;
      case BusToSignalAdapterPorts::ir:
      case BusToSignalAdapterPorts::ii:
        return SignalPortRole::Consumer;
      }
      return std::nullopt;
    }

    inline std::optional<SignalPortRole> portRole(Converter::RegcaPorts port)
    {
      switch (port)
      {
      case Converter::RegcaPorts::bus:
        return SignalPortRole::NonSignal;
      case Converter::RegcaPorts::ipcmd:
      case Converter::RegcaPorts::iqcmd:
        return SignalPortRole::Consumer;
      case Converter::RegcaPorts::ibranchr:
      case Converter::RegcaPorts::ibranchi:
      case Converter::RegcaPorts::pbranch:
      case Converter::RegcaPorts::qbranch:
        return SignalPortRole::Producer;
      }
      return std::nullopt;
    }

    inline std::optional<SignalPortRole> portRole(Converter::ReecbPorts port)
    {
      switch (port)
      {
      case Converter::ReecbPorts::bus:
        return SignalPortRole::NonSignal;
      case Converter::ReecbPorts::pe:
      case Converter::ReecbPorts::qgen:
      case Converter::ReecbPorts::qext:
      case Converter::ReecbPorts::pfaref:
      case Converter::ReecbPorts::pref:
        return SignalPortRole::Consumer;
      case Converter::ReecbPorts::iqcmd:
      case Converter::ReecbPorts::ipcmd:
        return SignalPortRole::Producer;
      }
      return std::nullopt;
    }

    inline std::optional<SignalPortRole> portRole(Converter::RepcaPorts port)
    {
      switch (port)
      {
      case Converter::RepcaPorts::bus:
        return SignalPortRole::NonSignal;
      case Converter::RepcaPorts::ibranchr:
      case Converter::RepcaPorts::ibranchi:
      case Converter::RepcaPorts::qbranch:
      case Converter::RepcaPorts::pbranch:
      case Converter::RepcaPorts::vref:
      case Converter::RepcaPorts::qref:
      case Converter::RepcaPorts::pplantref:
      case Converter::RepcaPorts::freq:
      case Converter::RepcaPorts::freqref:
        return SignalPortRole::Consumer;
      case Converter::RepcaPorts::qext:
      case Converter::RepcaPorts::pext:
        return SignalPortRole::Producer;
      }
      return std::nullopt;
    }

    inline std::optional<SignalPortRole> portRole(Exciter::Ieeet1Ports port)
    {
      switch (port)
      {
      case Exciter::Ieeet1Ports::bus:
        return SignalPortRole::NonSignal;
      case Exciter::Ieeet1Ports::speed:
      case Exciter::Ieeet1Ports::vs:
        return SignalPortRole::Consumer;
      case Exciter::Ieeet1Ports::efd:
        return SignalPortRole::Producer;
      }
      return std::nullopt;
    }

    inline std::optional<SignalPortRole> portRole(Exciter::SexsPtiPorts port)
    {
      switch (port)
      {
      case Exciter::SexsPtiPorts::bus:
        return SignalPortRole::NonSignal;
      case Exciter::SexsPtiPorts::vs:
        return SignalPortRole::Consumer;
      case Exciter::SexsPtiPorts::efd:
        return SignalPortRole::Producer;
      }
      return std::nullopt;
    }

    inline std::optional<SignalPortRole> portRole(Governor::Tgov1Ports port)
    {
      switch (port)
      {
      case Governor::Tgov1Ports::bus:
        return SignalPortRole::NonSignal;
      case Governor::Tgov1Ports::speed:
        return SignalPortRole::Consumer;
      case Governor::Tgov1Ports::pmech:
        return SignalPortRole::Producer;
      }
      return std::nullopt;
    }

    inline std::optional<SignalPortRole> portRole(Governor::HygovPorts port)
    {
      switch (port)
      {
      case Governor::HygovPorts::bus:
        return SignalPortRole::NonSignal;
      case Governor::HygovPorts::speed:
      case Governor::HygovPorts::pref:
      case Governor::HygovPorts::paux:
        return SignalPortRole::Consumer;
      case Governor::HygovPorts::pmech:
        return SignalPortRole::Producer;
      }
      return std::nullopt;
    }

    inline std::optional<SignalPortRole> portRole(Governor::GastPtiPorts port)
    {
      switch (port)
      {
      case Governor::GastPtiPorts::bus:
        return SignalPortRole::NonSignal;
      case Governor::GastPtiPorts::speed:
      case Governor::GastPtiPorts::pref:
        return SignalPortRole::Consumer;
      case Governor::GastPtiPorts::pmech:
        return SignalPortRole::Producer;
      }
      return std::nullopt;
    }

    inline std::optional<SignalPortRole> portRole(Stabilizer::IeeestPorts port)
    {
      switch (port)
      {
      case Stabilizer::IeeestPorts::input:
        return SignalPortRole::Consumer;
      case Stabilizer::IeeestPorts::output:
        return SignalPortRole::Producer;
      }
      return std::nullopt;
    }

    inline std::optional<SignalPortRole> portRole(GenrouPorts port)
    {
      switch (port)
      {
      case GenrouPorts::bus:
        return SignalPortRole::NonSignal;
      case GenrouPorts::pmech:
      case GenrouPorts::efd:
        return SignalPortRole::Consumer;
      case GenrouPorts::speed:
        return SignalPortRole::Producer;
      }
      return std::nullopt;
    }

    inline std::optional<SignalPortRole> portRole(GensalPorts port)
    {
      switch (port)
      {
      case GensalPorts::bus:
        return SignalPortRole::NonSignal;
      case GensalPorts::pmech:
      case GensalPorts::efd:
        return SignalPortRole::Consumer;
      case GensalPorts::speed:
        return SignalPortRole::Producer;
      }
      return std::nullopt;
    }

    struct ForcedOscillationInjection
    {
      std::string              id;
      std::string              target;
      std::string              mode{"add"};
      json                     params{json::object()};
      std::vector<std::string> mon{"in", "force", "out", "active"};
    };

    /**
     * @brief Describes an event that is used to modify the simulation at the
     * given time point
     */
    struct SystemEvent
    {
      /// Type of event determines action performed
      enum class Type
      {
        FAULT_ON,
        FAULT_OFF
      };

      /// Time event takes place
      double      time;
      /// Event type
      Type        type;
      /// ID of element used in event (e.g., bus fault id)
      std::size_t element_id;
    };

    /**
     * @brief Data defined in JSON file for parameterized study
     */
    struct StudyData
    {
      /// path to system model JSON file
      fs::path                                system_model_file;
      /// time step size
      double                                  dt;
      /// max time
      double                                  tmax;
      /// optional cap on IDA's internal step size
      std::optional<double>                   max_step;
      /// set of system events
      std::vector<SystemEvent>                events;
      /// forced oscillation injections for FO-only studies
      std::vector<ForcedOscillationInjection> forced_oscillations;
      /// path to output file
      fs::path                                output_file;
      /// path to reference file for validation
      fs::path                                reference_file;
      /// Error tolerance (between output file and reference file)
      double                                  error_tol;
      /// Instance of base model data
      SystemModelData<>                       model_data;
      /// Pristine model data, before any study-time injection transform
      SystemModelData<>                       base_model_data;
      /// True when the study uses forced_oscillations instead of events
      bool                                    forced_oscillation_study{false};

      bool isForcedOscillationStudy() const
      {
        return forced_oscillation_study;
      }
    };

    inline void from_json(const json& j, ForcedOscillationInjection& c)
    {
      j.at("id").get_to(c.id);
      j.at("target").get_to(c.target);
      c.mode = j.value("mode", std::string{"add"});

      if (j.contains("params"))
      {
        c.params = j.at("params");
        if (!c.params.is_object())
        {
          detail::studyError("Forced oscillation params must be an object");
        }
      }

      if (j.contains("mon"))
      {
        j.at("mon").get_to(c.mon);
      }
    }

    /**
     * @brief JSON parser implemntation for `StudyData`
     */
    inline void from_json(const json& j, StudyData& c)
    {
      using namespace magic_enum;

      static const std::set<std::string> known_keys = {
          "system_model_file",
          "dt",
          "tmax",
          "max_step",
          "output_file",
          "reference_file",
          "error_tolerance",
          "events",
          "forced_oscillations"};

      for (const auto& item : j.items())
      {
        const auto& key = item.key();
        if (!known_keys.contains(key))
        {
          detail::studyError("Unknown study field \"" + key + "\"");
        }
      }

      j.at("system_model_file").get_to(c.system_model_file);
      j.at("dt").get_to(c.dt);
      j.at("tmax").get_to(c.tmax);

      if (j.contains("max_step"))
      {
        c.max_step = j.at("max_step").get<double>();
        if (*c.max_step <= 0.0)
        {
          detail::studyError("max_step must be positive");
        }
      }

      if (j.contains("output_file"))
      {
        j.at("output_file").get_to(c.output_file);
      }

      if (j.contains("reference_file"))
      {
        j.at("reference_file").get_to(c.reference_file);
      }

      c.forced_oscillation_study = j.contains("forced_oscillations");
      if (c.forced_oscillation_study)
      {
        if (j.contains("events"))
        {
          detail::studyError("forced oscillation studies must not contain events");
        }
        j.at("forced_oscillations").get_to(c.forced_oscillations);
      }
      else
      {
        for (auto& raw_event : j.at("events"))
        {
          auto& event = c.events.emplace_back();
          raw_event.at("time").get_to(event.time);
          raw_event.at("element_id").get_to(event.element_id);

          auto type_str   = raw_event.at("type").get<std::string>();
          using EventType = SystemEvent::Type;
          auto type_wrap  = enum_cast<EventType>(type_str, case_insensitive);
          if (!type_wrap.has_value())
          {
            detail::studyError("Unable to parse event type \"" + type_str + "\"");
          }
          event.type = type_wrap.value();
        }
      }

      c.error_tol = j.value("error_tolerance", 1.0e-4);
    }

    namespace detail
    {
      struct ComponentTarget
      {
        std::string component_class;
        std::string component_id;
        std::string port;
      };

      struct ExplicitSignalTarget
      {
        std::size_t signal_id;
      };

      template <typename DataT>
      bool roleTableComplete()
      {
        for (auto port : magic_enum::enum_values<typename DataT::Ports>())
        {
          if (!portRole(port).has_value())
          {
            return false;
          }
        }
        return true;
      }

      template <typename DataT>
      DataT* findDevice(std::vector<DataT>& devices, const std::string& id)
      {
        auto iter = std::find_if(devices.begin(),
                                 devices.end(),
                                 [&id](const auto& device)
                                 { return device.disambiguation_string == id; });
        if (iter == devices.end())
        {
          return nullptr;
        }
        return &(*iter);
      }

      std::set<std::size_t> collectSignalIds(const SystemModelData<>& model)
      {
        std::set<std::size_t> ids;
        for (const auto& signal : model.signal)
        {
          ids.insert(signal.signal_id);
        }
        return ids;
      }

      std::size_t createSignal(
          SystemModelData<>&     model,
          std::set<std::size_t>& signal_ids,
          const std::string&     name)
      {
        std::size_t id = signal_ids.empty() ? 0 : (*signal_ids.rbegin() + 1);
        while (signal_ids.contains(id))
        {
          ++id;
        }

        SystemModelData<>::SignalDataT signal;
        signal.signal_id = id;
        signal.name      = name;
        model.signal.push_back(signal);
        signal_ids.insert(id);
        return id;
      }

      void ensureSignalExists(
          SystemModelData<>&     model,
          std::set<std::size_t>& signal_ids,
          std::size_t            signal_id,
          const std::string&     name)
      {
        if (signal_ids.contains(signal_id))
        {
          return;
        }

        SystemModelData<>::SignalDataT signal;
        signal.signal_id = signal_id;
        signal.name      = name;
        model.signal.push_back(signal);
        signal_ids.insert(signal_id);
      }

      template <typename DataT>
      bool vectorHasProducer(const std::vector<DataT>& devices, std::size_t signal_id)
      {
        for (const auto& device : devices)
        {
          for (const auto& [port, id] : device.ports)
          {
            auto role = portRole(port);
            if (role.has_value() && *role == SignalPortRole::Producer && id == signal_id)
            {
              return true;
            }
          }
        }
        return false;
      }

      bool signalHasProducer(const SystemModelData<>& model, std::size_t signal_id)
      {
        return vectorHasProducer(model.adapter, signal_id)
               || vectorHasProducer(model.regca, signal_id)
               || vectorHasProducer(model.reecb, signal_id)
               || vectorHasProducer(model.repca, signal_id)
               || vectorHasProducer(model.genrou, signal_id)
               || vectorHasProducer(model.gensal, signal_id)
               || vectorHasProducer(model.gov, signal_id)
               || vectorHasProducer(model.hygov, signal_id)
               || vectorHasProducer(model.gastpti, signal_id)
               || vectorHasProducer(model.exciter, signal_id)
               || vectorHasProducer(model.sexspti, signal_id)
               || vectorHasProducer(model.stabilizer, signal_id)
               || std::any_of(model.forced_oscillation.begin(),
                              model.forced_oscillation.end(),
                              [signal_id](const auto& forced)
                              {
                                using Ports = typename SystemModelData<>::ForcedOscillationDataT::Ports;
                                return forced.ports.contains(Ports::output)
                                       && forced.ports.at(Ports::output) == signal_id;
                              });
      }

      ComponentTarget parseComponentTarget(const std::string& target)
      {
        const auto colon = target.find(':');
        const auto dot   = target.rfind('.');
        if (colon == std::string::npos || dot == std::string::npos || dot <= colon + 1)
        {
          studyError("Malformed target \"" + target + "\"; expected Class:id.port or signal:ID");
        }

        ComponentTarget parsed;
        parsed.component_class = target.substr(0, colon);
        parsed.component_id    = target.substr(colon + 1, dot - colon - 1);
        parsed.port            = target.substr(dot + 1);
        if (parsed.component_class.empty() || parsed.component_id.empty() || parsed.port.empty())
        {
          studyError("Malformed target \"" + target + "\"; expected Class:id.port or signal:ID");
        }
        return parsed;
      }

      ExplicitSignalTarget parseSignalTarget(const std::string& target)
      {
        constexpr auto prefix = std::string_view{"signal:"};
        if (!target.starts_with(prefix))
        {
          studyError("Malformed signal target \"" + target + "\"; expected signal:ID");
        }

        const auto id_text = target.substr(prefix.size());
        if (id_text.empty())
        {
          studyError("Malformed signal target \"" + target + "\"; expected signal:ID");
        }

        try
        {
          return {static_cast<std::size_t>(std::stoull(id_text))};
        }
        catch (const std::exception&)
        {
          studyError("Malformed signal target \"" + target + "\"; expected numeric signal ID");
        }
      }

      void appendForcedOscillation(
          SystemModelData<>&                model,
          const ForcedOscillationInjection& injection,
          std::size_t                       output_signal,
          std::optional<std::size_t>        input_signal)
      {
        using DataT  = typename SystemModelData<>::ForcedOscillationDataT;
        using Params = typename DataT::Parameters;
        using Ports  = typename DataT::Ports;
        using Mon    = typename DataT::MonitorableVariables;

        DataT data;
        data.device_class          = "ForcedOscillation";
        data.disambiguation_string = injection.id;
        data.ports[Ports::output]  = output_signal;
        if (input_signal.has_value())
        {
          data.ports[Ports::input] = *input_signal;
        }

        for (const auto& [name, value] : injection.params.items())
        {
          auto parameter = magic_enum::enum_cast<Params>(name);
          if (!parameter.has_value())
          {
            studyError("Unknown ForcedOscillation parameter \"" + name + "\" in injection \"" + injection.id + "\"");
          }
          if (!value.is_number())
          {
            studyError("ForcedOscillation parameter \"" + name + "\" in injection \"" + injection.id + "\" must be numeric");
          }
          data.parameters[*parameter] = value.get<double>();
        }

        for (const auto& variable : injection.mon)
        {
          auto monitored = magic_enum::enum_cast<Mon>(variable);
          if (!monitored.has_value())
          {
            studyError("Unknown ForcedOscillation monitor \"" + variable + "\" in injection \"" + injection.id + "\"");
          }
          data.monitored_variables.insert(*monitored);
        }

        model.forced_oscillation.push_back(data);
      }

      template <typename DataT>
      void applyComponentInjection(
          SystemModelData<>&                model,
          std::vector<DataT>&               devices,
          const ComponentTarget&            target,
          const ForcedOscillationInjection& injection,
          std::set<std::size_t>&            signal_ids)
      {
        using Ports = typename DataT::Ports;

        auto port = magic_enum::enum_cast<Ports>(target.port);
        if (!port.has_value())
        {
          studyError("Unknown port \"" + target.port + "\" on target \"" + injection.target + "\"");
        }

        auto role = portRole(*port);
        if (!role.has_value())
        {
          studyError("No signal-port role is defined for target \"" + injection.target + "\"");
        }

        if (*role == SignalPortRole::Producer)
        {
          studyError("Target \"" + injection.target + "\" is a producer port; target the downstream consumer port instead");
        }
        if (*role == SignalPortRole::NonSignal)
        {
          studyError("Target \"" + injection.target + "\" is not a signal consumer port");
        }

        auto* device = findDevice(devices, target.component_id);
        if (device == nullptr)
        {
          studyError("Unable to find " + target.component_class + " device \"" + target.component_id + "\"");
        }

        const bool is_add   = injection.mode == "add";
        const bool is_drive = injection.mode == "drive";
        if (!is_add && !is_drive)
        {
          studyError("Unknown forced oscillation mode \"" + injection.mode + "\"");
        }

        std::optional<std::size_t> input_signal;
        if (is_add)
        {
          if (!device->ports.contains(*port))
          {
            studyError("add mode requires an existing consumer signal on target \"" + injection.target + "\"");
          }
          input_signal = device->ports.at(*port);
        }

        const auto output_signal = createSignal(model, signal_ids, injection.id + "_out");
        device->ports[*port]     = output_signal;
        appendForcedOscillation(model, injection, output_signal, input_signal);
      }

      void applyComponentInjection(
          SystemModelData<>&                model,
          const ComponentTarget&            target,
          const ForcedOscillationInjection& injection,
          std::set<std::size_t>&            signal_ids)
      {
        if (target.component_class == "BusToSignalAdapter")
        {
          applyComponentInjection(model, model.adapter, target, injection, signal_ids);
        }
        else if (target.component_class == "Regca")
        {
          applyComponentInjection(model, model.regca, target, injection, signal_ids);
        }
        else if (target.component_class == "Reecb")
        {
          applyComponentInjection(model, model.reecb, target, injection, signal_ids);
        }
        else if (target.component_class == "Repca")
        {
          applyComponentInjection(model, model.repca, target, injection, signal_ids);
        }
        else if (target.component_class == "Ieeet1")
        {
          applyComponentInjection(model, model.exciter, target, injection, signal_ids);
        }
        else if (target.component_class == "SexsPti")
        {
          applyComponentInjection(model, model.sexspti, target, injection, signal_ids);
        }
        else if (target.component_class == "Tgov1")
        {
          applyComponentInjection(model, model.gov, target, injection, signal_ids);
        }
        else if (target.component_class == "Hygov")
        {
          applyComponentInjection(model, model.hygov, target, injection, signal_ids);
        }
        else if (target.component_class == "GastPti")
        {
          applyComponentInjection(model, model.gastpti, target, injection, signal_ids);
        }
        else if (target.component_class == "Ieeest")
        {
          applyComponentInjection(model, model.stabilizer, target, injection, signal_ids);
        }
        else if (target.component_class == "Genrou")
        {
          applyComponentInjection(model, model.genrou, target, injection, signal_ids);
        }
        else if (target.component_class == "Gensal")
        {
          applyComponentInjection(model, model.gensal, target, injection, signal_ids);
        }
        else
        {
          studyError("Unknown forced oscillation target class \"" + target.component_class + "\"");
        }
      }

      double readInjectionParameter(
          const ForcedOscillationInjection& injection,
          const std::string&                name,
          double                            fallback)
      {
        if (!injection.params.contains(name))
        {
          return fallback;
        }
        const auto& value = injection.params.at(name);
        if (!value.is_number())
        {
          studyError("ForcedOscillation parameter \"" + name + "\" in injection \"" + injection.id + "\" must be numeric");
        }
        return value.get<double>();
      }
    } // namespace detail

    bool forcedOscillationPortRoleTableComplete()
    {
      return detail::roleTableComplete<SystemModelData<>::BusToSignalAdapterDataT>()
             && detail::roleTableComplete<SystemModelData<>::RegcaDataT>()
             && detail::roleTableComplete<SystemModelData<>::ReecbDataT>()
             && detail::roleTableComplete<SystemModelData<>::RepcaDataT>()
             && detail::roleTableComplete<SystemModelData<>::Ieeet1DataT>()
             && detail::roleTableComplete<SystemModelData<>::SexsPtiDataT>()
             && detail::roleTableComplete<SystemModelData<>::Tgov1DataT>()
             && detail::roleTableComplete<SystemModelData<>::HygovDataT>()
             && detail::roleTableComplete<SystemModelData<>::GastPtiDataT>()
             && detail::roleTableComplete<SystemModelData<>::IeeestDataT>()
             && detail::roleTableComplete<SystemModelData<>::GenrouDataT>()
             && detail::roleTableComplete<SystemModelData<>::GensalDataT>();
    }

    SystemModelData<> applyInjections(
        const SystemModelData<>&                       base,
        const std::vector<ForcedOscillationInjection>& injections)
    {
      auto model      = base;
      auto signal_ids = detail::collectSignalIds(model);

      std::set<std::string> ids;
      for (const auto& forced : model.forced_oscillation)
      {
        ids.insert(forced.disambiguation_string);
      }

      for (const auto& injection : injections)
      {
        if (!ids.insert(injection.id).second)
        {
          detail::studyError("Duplicate forced oscillation injection id \"" + injection.id + "\"");
        }

        if (injection.mode == "signal")
        {
          auto target = detail::parseSignalTarget(injection.target);
          if (detail::signalHasProducer(model, target.signal_id))
          {
            detail::studyError("signal target \"" + injection.target + "\" already has a producer");
          }
          detail::ensureSignalExists(model, signal_ids, target.signal_id, injection.id + "_out");
          detail::appendForcedOscillation(model, injection, target.signal_id, std::nullopt);
        }
        else
        {
          auto target = detail::parseComponentTarget(injection.target);
          detail::applyComponentInjection(model, target, injection, signal_ids);
        }
      }

      return model;
    }

    double forcedOscillationMaxFrequency(
        const std::vector<ForcedOscillationInjection>& injections,
        double                                         tmax)
    {
      double max_frequency = 0.0;
      for (const auto& injection : injections)
      {
        const double f   = detail::readInjectionParameter(injection, "f", 0.0);
        const double Kf  = detail::readInjectionParameter(injection, "Kf", 0.0);
        const double Ton = detail::readInjectionParameter(injection, "Ton", 0.0);
        const double tau = std::max(tmax - Ton, 0.0);
        max_frequency    = std::max(max_frequency, f + Kf * tau);
      }
      return std::max(max_frequency, 0.0);
    }

    double defaultForcedOscillationMaxStep(
        const std::vector<ForcedOscillationInjection>& injections,
        double                                         dt,
        double                                         tmax)
    {
      const double f_max = forcedOscillationMaxFrequency(injections, tmax);
      if (f_max == 0.0)
      {
        return dt;
      }
      return std::min(dt, 1.0 / (20.0 * f_max));
    }

    /**
     * @brief Check for existence and successful input file open
     */
    std::ifstream openFile(const fs::path& file_path)
    {
      if (!exists(file_path))
      {
        detail::studyError("File not found: " + file_path.string());
      }
      auto fs = std::ifstream(file_path);
      if (!fs)
      {
        detail::studyError("Failed to open file: " + file_path.string());
      }
      return fs;
    }

    /**
     * @brief Wrapper function to parse `StudyData` from JSON and perform
     * follow-up configuration
     */
    StudyData parseStudyData(const fs::path& file_path)
    {
      auto data = StudyData(json::parse(openFile(file_path)));

      auto loc = file_path.parent_path();
      if (!data.system_model_file.is_absolute())
      {
        data.system_model_file = loc / data.system_model_file;
      }
      if (!data.reference_file.empty())
      {
        if (!data.reference_file.is_absolute())
        {
          data.reference_file = loc / data.reference_file;
        }
      }

      auto csv        = ::GridKit::Model::VariableMonitorFormat::CSV;
      data.model_data = parseSystemModelData(data.system_model_file);
      std::string model_output_file;
      // Find output file (CSV) specified in model input file
      for (const auto& sink : data.model_data.monitor_sink)
      {
        if (sink.format == csv && sink.delim == ",")
        {
          model_output_file = sink.file_name;
        }
      }

      if (model_output_file.empty())
      {
        // Add study output file to model if one did not already exist
        data.model_data.monitor_sink.emplace_back(csv, data.output_file);
      }
      else
      {
        if (data.output_file.empty())
        {
          data.output_file = model_output_file;
        }
        else
        {
          // If model file already specifies a CSV output file, then the study
          // output file must be a symlink to the model output file
          if (exists(data.output_file))
          {
            if ((!is_symlink(data.output_file)) || (read_symlink(data.output_file) != model_output_file))
            {
              detail::studyError("Study output file not usable");
            }
          }
          else
          {
            fs::create_symlink(model_output_file, data.output_file);
          }
        }
      }

      data.base_model_data = data.model_data;

      if (data.isForcedOscillationStudy() && !data.max_step.has_value())
      {
        data.max_step = defaultForcedOscillationMaxStep(data.forced_oscillations, data.dt, data.tmax);
        if (*data.max_step < data.dt)
        {
          Log::warning() << "Forced oscillation max_step set to " << *data.max_step
                         << " to resolve f_max = "
                         << forcedOscillationMaxFrequency(data.forced_oscillations, data.tmax)
                         << '\n';
        }
      }

      return data;
    }

    void checkCommandLine(int argc, const std::string& appName)
    {
      if (argc < 2)
      {
        Log::error() << "No input file provided" << std::endl;
        std::cout << std::format(
            "\n"
            "Usage:\n"
            "       {} <json-input-file>\n"
            "\n"
            "Please provide a json input file for the study to run.\n"
            "\n",
            appName);
        exit(1);
      }
    }

    Testing::TestStatus checkErrors(
        const StudyData& study_data,
        bool             print_results = true)
    {
      // Generate aggregate errors comparing variable output to reference solution
      auto func   = std::string{"monitor file vs reference file"};
      auto status = Testing::TestStatus{func.c_str()};

      const auto& out_file = study_data.output_file;
      const auto& ref_file = study_data.reference_file;
      if (!out_file.empty() && !ref_file.empty())
      {
        auto errorSet = Testing::compareCSV(out_file, ref_file);

        // Print the errors
        if (print_results)
        {
          errorSet.display();
        }

        // Check against specified tolerance
        status *= errorSet.total.max_error < study_data.error_tol;

        if (print_results)
        {
          status.report();
        }
      }
      return status;
    }

  } // namespace PhasorDynamics
} // namespace GridKit
