/**
 * @file StateData.cpp
 * @brief JSON input and output for `StateData`.
 */

#include <fstream>
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

      StateRecord readRecord(const std::string& id, const json& object)
      {
        StateRecord record;
        for (const auto& [key, value] : object.items())
        {
          if (value.is_boolean())
          {
            record.flags[key] = value.get<bool>();
          }
          else if (value.is_number())
          {
            record.values[key] = value.get<double>();
          }
          else if (!value.is_null())
          {
            throw std::invalid_argument("State entry " + id + "." + key + " must be a number or a boolean");
          }
        }
        return record;
      }

      std::map<std::string, StateRecord> readRecords(const json& document, const char* section)
      {
        std::map<std::string, StateRecord> records;
        if (document.contains(section))
        {
          for (const auto& [id, object] : document.at(section).items())
          {
            records.emplace(id, readRecord(id, object));
          }
        }
        return records;
      }

      json writeRecords(const std::map<std::string, StateRecord>& records)
      {
        json section = json::object();
        for (const auto& [id, record] : records)
        {
          json object = json::object();
          for (const auto& [key, value] : record.values)
          {
            object[key] = value;
          }
          for (const auto& [key, flag] : record.flags)
          {
            object[key] = flag;
          }
          section[id] = object;
        }
        return section;
      }

      template <typename ValueT>
      void readOptional(const json& object, const char* key, std::optional<ValueT>& value)
      {
        if (object.contains(key) && !object.at(key).is_null())
        {
          value = object.at(key).get<ValueT>();
        }
      }

      template <typename ValueT>
      void writeOptional(json& object, const char* key, const std::optional<ValueT>& value)
      {
        if (value)
        {
          object[key] = *value;
        }
      }
    } // namespace

    double StateRecord::value(const std::string& key, double fallback) const
    {
      const auto entry = values.find(key);
      if (entry == values.end())
      {
        return fallback;
      }
      return entry->second;
    }

    bool StateRecord::flag(const std::string& key, bool fallback) const
    {
      const auto entry = flags.find(key);
      if (entry == flags.end())
      {
        return fallback;
      }
      return entry->second;
    }

    const StateRecord* StateData::bus(std::size_t number) const
    {
      const auto entry = buses.find(busKey(number));
      if (entry == buses.end())
      {
        return nullptr;
      }
      return &entry->second;
    }

    const StateRecord* StateData::device(const std::string& id) const
    {
      const auto entry = devices.find(id);
      if (entry == devices.end())
      {
        return nullptr;
      }
      return &entry->second;
    }

    std::string busKey(std::size_t number)
    {
      return "bus_id_" + std::to_string(number);
    }

    std::pair<std::string, std::string> currentKeys(std::size_t terminal, std::size_t terminal_count)
    {
      if (terminal_count == 1)
      {
        return {"ir", "ii"};
      }
      const auto suffix = std::to_string(terminal + 1);
      return {"ir" + suffix, "ii" + suffix};
    }

    bool terminalPower(const StateData&   state,
                       const std::string& id,
                       std::size_t        number,
                       std::size_t        terminal,
                       std::size_t        terminal_count,
                       double&            p,
                       double&            q)
    {
      const StateRecord* device = state.device(id);
      const StateRecord* bus    = state.bus(number);
      if (device == nullptr || bus == nullptr)
      {
        return false;
      }

      const auto [ir_key, ii_key] = currentKeys(terminal, terminal_count);
      if (!device->values.contains(ir_key) || !device->values.contains(ii_key)
          || !bus->values.contains("vr") || !bus->values.contains("vi"))
      {
        return false;
      }

      const double vr = bus->values.at("vr");
      const double vi = bus->values.at("vi");
      const double ir = device->values.at(ir_key);
      const double ii = device->values.at(ii_key);

      p = vr * ir + vi * ii;
      q = vi * ir - vr * ii;
      return true;
    }

    void setTerminalCurrent(StateData&         state,
                            const std::string& id,
                            std::size_t        number,
                            std::size_t        terminal,
                            std::size_t        terminal_count,
                            double             p,
                            double             q)
    {
      const auto&  voltage = state.buses.at(busKey(number)).values;
      const double vr      = voltage.at("vr");
      const double vi      = voltage.at("vi");
      const double v2      = vr * vr + vi * vi;

      // I = conj((P + jQ) / V)
      const auto [ir_key, ii_key] = currentKeys(terminal, terminal_count);
      auto& current               = state.devices[id].values;
      current[ir_key]             = (p * vr + q * vi) / v2;
      current[ii_key]             = (p * vi - q * vr) / v2;
    }

    StateData parseStateData(std::istream& stream)
    {
      const json document = json::parse(stream);

      StateData state;
      if (document.contains("header"))
      {
        const auto& header = document.at("header");
        readOptional(header, "version", state.header.version);
        readOptional(header, "time", state.header.time);
        readOptional(header, "created", state.header.created);
        readOptional(header, "description", state.header.description);
      }
      state.buses   = readRecords(document, "buses");
      state.devices = readRecords(document, "devices");
      return state;
    }

    StateData parseStateData(const std::filesystem::path& file)
    {
      std::ifstream stream(file);
      if (!stream)
      {
        throw std::runtime_error("Could not open state file " + file.string());
      }
      return parseStateData(stream);
    }

    void writeStateData(const StateData& state, std::ostream& stream)
    {
      json header = json::object();
      writeOptional(header, "version", state.header.version);
      writeOptional(header, "time", state.header.time);
      writeOptional(header, "created", state.header.created);
      writeOptional(header, "description", state.header.description);

      json document       = json::object();
      document["header"]  = header;
      document["buses"]   = writeRecords(state.buses);
      document["devices"] = writeRecords(state.devices);

      stream << document.dump(4) << "\n";
    }

    void writeStateData(const StateData& state, const std::filesystem::path& file)
    {
      std::ofstream stream(file);
      if (!stream)
      {
        throw std::runtime_error("Could not open state file " + file.string());
      }
      writeStateData(state, stream);
    }
  } // namespace Model
} // namespace GridKit
