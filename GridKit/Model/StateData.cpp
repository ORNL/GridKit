#include <fstream>
#include <sstream>
#include <stdexcept>

#include <nlohmann/json.hpp>

#include <GridKit/Model/StateData.hpp>

namespace GridKit
{
  namespace Model
  {
    namespace
    {
      using json = nlohmann::json;

      template <typename ValueT>
      void parseOptional(const json& object, const char* key, std::optional<ValueT>& value)
      {
        auto entry = object.find(key);
        if (entry != object.end() && !entry->is_null())
        {
          value = entry->get<ValueT>();
        }
      }

      void requireObject(const json& value, const char* name)
      {
        if (!value.is_object())
        {
          std::stringstream message;
          message << name << " must be a JSON object";
          throw std::invalid_argument(message.str());
        }
      }

      InjectionState parseInjectionState(const json& value)
      {
        requireObject(value, "Injection state");

        InjectionState state;
        parseOptional(value, "ir", state.ir);
        parseOptional(value, "ii", state.ii);
        return state;
      }

      BusState parseBusState(const json& value)
      {
        requireObject(value, "Bus state");

        BusState state;
        parseOptional(value, "vr", state.vr);
        parseOptional(value, "vi", state.vi);

        auto injections = value.find("injections");
        if (injections != value.end() && !injections->is_null())
        {
          requireObject(*injections, "Bus injections");
          for (const auto& [id, injection] : injections->items())
          {
            if (!injection.is_null())
            {
              state.injections.emplace(id, parseInjectionState(injection));
            }
          }
        }
        return state;
      }

      DeviceState parseDeviceState(const json& value)
      {
        requireObject(value, "Device state");

        DeviceState state;
        parseOptional(value, "online", state.online);
        parseOptional(value, "open", state.open);
        parseOptional(value, "p", state.p);
        parseOptional(value, "q", state.q);
        parseOptional(value, "tap", state.tap);
        parseOptional(value, "phase", state.phase);
        return state;
      }

      StateData parseStateData(const json& value)
      {
        requireObject(value, "State data");

        StateData state;

        auto header = value.find("header");
        if (header != value.end() && !header->is_null())
        {
          requireObject(*header, "State header");
          state.header.emplace();
          parseOptional(*header, "version", state.header->version);
          parseOptional(*header, "time", state.header->time);
          parseOptional(*header, "created", state.header->created);
          parseOptional(*header, "description", state.header->description);
        }

        auto buses = value.find("buses");
        if (buses != value.end() && !buses->is_null())
        {
          requireObject(*buses, "State buses");
          for (const auto& [id, bus] : buses->items())
          {
            if (!bus.is_null())
            {
              state.buses.emplace(id, parseBusState(bus));
            }
          }
        }

        auto devices = value.find("devices");
        if (devices != value.end() && !devices->is_null())
        {
          requireObject(*devices, "State devices");
          for (const auto& [id, device] : devices->items())
          {
            if (!device.is_null())
            {
              state.devices.emplace(id, parseDeviceState(device));
            }
          }
        }

        return state;
      }
    } // namespace

    StateData parseStateData(std::istream& stream)
    {
      return parseStateData(json::parse(stream));
    }

    StateData parseStateData(const std::filesystem::path& file_path)
    {
      std::ifstream stream(file_path);
      if (!stream)
      {
        std::stringstream message;
        message << "Could not open file: " << file_path;
        throw std::runtime_error(message.str());
      }
      return parseStateData(stream);
    }
  } // namespace Model
} // namespace GridKit
