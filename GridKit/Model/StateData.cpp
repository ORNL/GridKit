#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>

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
        parseOptional(value, "ia", state.ia);
        parseOptional(value, "ib", state.ib);
        parseOptional(value, "ic", state.ic);
        return state;
      }

      BusState parseBusState(const json& value)
      {
        requireObject(value, "Bus state");

        BusState state;
        parseOptional(value, "vr", state.vr);
        parseOptional(value, "vi", state.vi);
        parseOptional(value, "va", state.va);
        parseOptional(value, "vb", state.vb);
        parseOptional(value, "vc", state.vc);

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
        parseOptional(value, "active", state.active);
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

      template <typename ValueT>
      void writeOptional(json&                        object,
                         const char*                  key,
                         const std::optional<ValueT>& value)
      {
        if (value.has_value())
        {
          if constexpr (std::is_floating_point_v<ValueT>)
          {
            if (!std::isfinite(*value))
            {
              std::stringstream message;
              message << "State field \"" << key << "\" must be finite";
              throw std::invalid_argument(message.str());
            }
          }
          object[key] = *value;
        }
      }

      json writeInjectionState(const InjectionState& state)
      {
        json value = json::object();
        writeOptional(value, "ir", state.ir);
        writeOptional(value, "ii", state.ii);
        writeOptional(value, "ia", state.ia);
        writeOptional(value, "ib", state.ib);
        writeOptional(value, "ic", state.ic);
        return value;
      }

      json writeBusState(const BusState& state)
      {
        json value = json::object();
        writeOptional(value, "vr", state.vr);
        writeOptional(value, "vi", state.vi);
        writeOptional(value, "va", state.va);
        writeOptional(value, "vb", state.vb);
        writeOptional(value, "vc", state.vc);

        json injections = json::object();
        for (const auto& [id, injection] : state.injections)
        {
          injections[id] = writeInjectionState(injection);
        }
        value["injections"] = std::move(injections);
        return value;
      }

      json writeDeviceState(const DeviceState& state)
      {
        json value = json::object();
        writeOptional(value, "active", state.active);
        writeOptional(value, "online", state.online);
        writeOptional(value, "open", state.open);
        writeOptional(value, "p", state.p);
        writeOptional(value, "q", state.q);
        writeOptional(value, "tap", state.tap);
        writeOptional(value, "phase", state.phase);
        return value;
      }

      json serializeStateData(const StateData& state)
      {
        json value = json::object();

        if (state.header.has_value())
        {
          json header = json::object();
          writeOptional(header, "version", state.header->version);
          writeOptional(header, "time", state.header->time);
          writeOptional(header, "created", state.header->created);
          writeOptional(header, "description", state.header->description);
          value["header"] = std::move(header);
        }

        json buses = json::object();
        for (const auto& [id, bus] : state.buses)
        {
          buses[id] = writeBusState(bus);
        }
        value["buses"] = std::move(buses);

        json devices = json::object();
        for (const auto& [id, device] : state.devices)
        {
          devices[id] = writeDeviceState(device);
        }
        value["devices"] = std::move(devices);

        return value;
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

    void writeStateData(const StateData& state, std::ostream& stream)
    {
      const json value = serializeStateData(state);
      stream << value.dump(4) << '\n';
      if (!stream)
      {
        throw std::runtime_error("Could not write state data");
      }
    }

    void writeStateData(const StateData&             state,
                        const std::filesystem::path& file_path)
    {
      const json    value = serializeStateData(state);
      std::ofstream stream(file_path);
      if (!stream)
      {
        std::stringstream message;
        message << "Could not open file: " << file_path;
        throw std::runtime_error(message.str());
      }
      stream << value.dump(4) << '\n';
      if (!stream)
      {
        throw std::runtime_error("Could not write state data");
      }
    }
  } // namespace Model
} // namespace GridKit
