#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>

#include <nlohmann/json.hpp>

#include <GridKit/Model/OPF/SystemData.hpp>

namespace GridKit
{
  namespace OPF
  {
    namespace
    {
      using Json           = nlohmann::json;
      using SystemDataT    = SystemData<>;
      using RealT          = SystemDataT::RealT;
      using IdxT           = SystemDataT::IdxT;
      using BusDataT       = SystemDataT::BusDataT;
      using BranchDataT    = SystemDataT::BranchDataT;
      using GeneratorDataT = SystemDataT::GeneratorDataT;
      using LoadDataT      = SystemDataT::LoadDataT;
      using ShuntDataT     = SystemDataT::ShuntDataT;
      using DeviceDataT    = SystemDataT::DeviceDataT;

      constexpr unsigned int SUPPORTED_FORMAT_VERSION  = 0;
      constexpr unsigned int SUPPORTED_FORMAT_REVISION = 1;

      [[noreturn]] void throwParseError(const std::string& context,
                                        const std::string& message)
      {
        throw std::invalid_argument(context + ": " + message);
      }

      void requireObject(const Json& value, const std::string& context)
      {
        if (!value.is_object())
        {
          throwParseError(context,
                          "expected an object; got " + std::string(value.type_name()));
        }
      }

      void requireArray(const Json& value, const std::string& context)
      {
        if (!value.is_array())
        {
          throwParseError(context,
                          "expected an array; got " + std::string(value.type_name()));
        }
      }

      bool isKnownField(std::string_view                        field,
                        std::initializer_list<std::string_view> known_fields)
      {
        for (const auto known_field : known_fields)
        {
          if (field == known_field)
          {
            return true;
          }
        }
        return false;
      }

      void rejectUnknownFields(const Json&                             object,
                               std::initializer_list<std::string_view> known_fields,
                               const std::string&                      context)
      {
        requireObject(object, context);
        for (const auto& entry : object.items())
        {
          const auto& field = entry.key();
          if (!isKnownField(field, known_fields))
          {
            throwParseError(context + "." + field, "unknown field");
          }
        }
      }

      const Json& requireMember(const Json&        object,
                                const char*        field,
                                const std::string& context)
      {
        requireObject(object, context);
        const auto entry = object.find(field);
        if (entry == object.end())
        {
          throwParseError(context + "." + field, "required field is missing");
        }
        if (entry->is_null())
        {
          throwParseError(context + "." + field, "required field must not be null");
        }
        return *entry;
      }

      const Json& requireObjectMember(const Json&        object,
                                      const char*        field,
                                      const std::string& context)
      {
        const Json& value = requireMember(object, field, context);
        requireObject(value, context + "." + field);
        return value;
      }

      const Json& requireArrayMember(const Json&        object,
                                     const char*        field,
                                     const std::string& context)
      {
        const Json& value = requireMember(object, field, context);
        requireArray(value, context + "." + field);
        return value;
      }

      std::string requireString(const Json&        object,
                                const char*        field,
                                const std::string& context)
      {
        const Json& value         = requireMember(object, field, context);
        const auto  field_context = context + "." + field;
        if (!value.is_string())
        {
          throwParseError(field_context,
                          "expected a string; got " + std::string(value.type_name()));
        }
        return value.get<std::string>();
      }

      RealT parseReal(const Json& value, const std::string& context)
      {
        if (!value.is_number())
        {
          throwParseError(context,
                          "expected a number; got " + std::string(value.type_name()));
        }

        RealT result{};
        try
        {
          result = value.get<RealT>();
        }
        catch (const Json::exception& error)
        {
          throwParseError(context, "invalid number: " + std::string(error.what()));
        }

        if (!std::isfinite(result))
        {
          throwParseError(context, "expected a finite number");
        }
        return result;
      }

      RealT requireReal(const Json&        object,
                        const char*        field,
                        const std::string& context)
      {
        return parseReal(requireMember(object, field, context), context + "." + field);
      }

      void readOptionalReal(const Json&           object,
                            const char*           field,
                            std::optional<RealT>& result,
                            const std::string&    context)
      {
        const auto entry = object.find(field);
        if (entry == object.end() || entry->is_null())
        {
          result.reset();
          return;
        }
        result = parseReal(*entry, context + "." + field);
      }

      RealT readRealOrDefault(const Json&        object,
                              const char*        field,
                              RealT              default_value,
                              const std::string& context)
      {
        const auto entry = object.find(field);
        if (entry == object.end())
        {
          return default_value;
        }
        if (entry->is_null())
        {
          throwParseError(context + "." + field, "field must not be null");
        }
        return parseReal(*entry, context + "." + field);
      }

      void readOptionalString(const Json&                 object,
                              const char*                 field,
                              std::optional<std::string>& result,
                              const std::string&          context)
      {
        const auto entry = object.find(field);
        if (entry == object.end() || entry->is_null())
        {
          result.reset();
          return;
        }
        if (!entry->is_string())
        {
          throwParseError(context + "." + field,
                          "expected a string or null; got "
                              + std::string(entry->type_name()));
        }
        result = entry->get<std::string>();
      }

      template <typename UnsignedT>
      UnsignedT parseUnsignedInteger(const Json& value, const std::string& context)
      {
        static_assert(std::is_integral_v<UnsignedT> && std::is_unsigned_v<UnsignedT>);

        std::uint64_t result{};
        try
        {
          if (value.is_number_unsigned())
          {
            result = value.get<std::uint64_t>();
          }
          else if (value.is_number_integer())
          {
            const auto signed_value = value.get<std::int64_t>();
            if (signed_value < 0)
            {
              throwParseError(context, "expected a non-negative integer");
            }
            result = static_cast<std::uint64_t>(signed_value);
          }
          else
          {
            throwParseError(context,
                            "expected a non-negative integer; got "
                                + std::string(value.type_name()));
          }
        }
        catch (const Json::exception& error)
        {
          throwParseError(context, "invalid integer: " + std::string(error.what()));
        }

        if (result > static_cast<std::uint64_t>(std::numeric_limits<UnsignedT>::max()))
        {
          throwParseError(context, "integer is outside the supported range");
        }
        return static_cast<UnsignedT>(result);
      }

      template <typename UnsignedT>
      UnsignedT requireUnsignedInteger(const Json&        object,
                                       const char*        field,
                                       const std::string& context)
      {
        return parseUnsignedInteger<UnsignedT>(requireMember(object, field, context),
                                               context + "." + field);
      }

      CaseHeader parseHeader(const Json& value)
      {
        const std::string context = "header";
        rejectUnknownFields(value,
                            {"format_version",
                             "format_revision",
                             "case_name",
                             "case_date_time",
                             "case_description",
                             "case_comments"},
                            context);

        CaseHeader header;
        header.format_version = requireUnsignedInteger<unsigned int>(
            value, "format_version", context);
        header.format_revision = requireUnsignedInteger<unsigned int>(
            value, "format_revision", context);
        if (header.format_version != SUPPORTED_FORMAT_VERSION)
        {
          throwParseError(context + ".format_version",
                          "unsupported format version "
                              + std::to_string(header.format_version)
                              + "; expected "
                              + std::to_string(SUPPORTED_FORMAT_VERSION));
        }
        if (header.format_revision != SUPPORTED_FORMAT_REVISION)
        {
          throwParseError(context + ".format_revision",
                          "unsupported format revision "
                              + std::to_string(header.format_revision)
                              + "; expected "
                              + std::to_string(SUPPORTED_FORMAT_REVISION));
        }
        header.case_name = requireString(value, "case_name", context);
        readOptionalString(value, "case_date_time", header.case_date_time, context);
        readOptionalString(value, "case_description", header.case_description, context);
        readOptionalString(value, "case_comments", header.case_comments, context);
        return header;
      }

      SystemParameters<RealT> parseSystemParameters(const Json& value)
      {
        const std::string context = "params";
        rejectUnknownFields(value, {"freq_base", "va_base"}, context);

        SystemParameters<RealT> params;
        params.freq_base = requireReal(value, "freq_base", context);
        params.va_base   = requireReal(value, "va_base", context);
        return params;
      }

      BusDataT parseBus(const Json& value, std::size_t position)
      {
        const auto base_context = "buses[" + std::to_string(position) + "]";
        rejectUnknownFields(value, {"class", "id", "params"}, base_context);

        const auto  bus_class = requireString(value, "class", base_context);
        const auto  id        = requireString(value, "id", base_context);
        const auto  context   = base_context + " (" + bus_class + " \"" + id + "\")";
        const auto& params    = requireObjectMember(value, "params", context);
        rejectUnknownFields(params,
                            {"number", "kv", "vmin", "vmax"},
                            context + ".params");

        BusDataT bus;
        if (bus_class == "Bus")
        {
          bus.bus_class = BusClass::BUS;
        }
        else if (bus_class == "Slack")
        {
          bus.bus_class = BusClass::SLACK;
        }
        else
        {
          throwParseError(base_context + ".class",
                          "unsupported bus class \"" + bus_class
                              + "\"; expected \"Bus\" or \"Slack\"");
        }

        bus.id     = id;
        bus.number = requireUnsignedInteger<IdxT>(params, "number", context + ".params");
        bus.kv     = requireReal(params, "kv", context + ".params");
        readOptionalReal(params, "vmin", bus.vmin, context + ".params");
        readOptionalReal(params, "vmax", bus.vmax, context + ".params");
        return bus;
      }

      BranchDataT parseBranch(const Json&        buses,
                              const Json&        params,
                              const std::string& id,
                              const std::string& context)
      {
        rejectUnknownFields(buses, {"from", "to"}, context + ".buses");
        rejectUnknownFields(params,
                            {"R", "X", "G", "B", "smax"},
                            context + ".params");

        BranchDataT branch;
        branch.id   = id;
        branch.from = requireUnsignedInteger<IdxT>(buses, "from", context + ".buses");
        branch.to   = requireUnsignedInteger<IdxT>(buses, "to", context + ".buses");
        branch.R    = requireReal(params, "R", context + ".params");
        branch.X    = requireReal(params, "X", context + ".params");
        branch.G    = requireReal(params, "G", context + ".params");
        branch.B    = requireReal(params, "B", context + ".params");
        readOptionalReal(params, "smax", branch.smax, context + ".params");
        return branch;
      }

      GeneratorDataT parseGenerator(const Json&        buses,
                                    const Json&        params,
                                    const std::string& id,
                                    const std::string& context)
      {
        rejectUnknownFields(buses, {"bus"}, context + ".buses");
        rejectUnknownFields(
            params,
            {"pmin", "pmax", "qmin", "qmax", "mva", "c0", "c1", "c2"},
            context + ".params");

        GeneratorDataT generator;
        generator.id  = id;
        generator.bus = requireUnsignedInteger<IdxT>(buses, "bus", context + ".buses");
        readOptionalReal(params, "pmin", generator.pmin, context + ".params");
        readOptionalReal(params, "pmax", generator.pmax, context + ".params");
        readOptionalReal(params, "qmin", generator.qmin, context + ".params");
        readOptionalReal(params, "qmax", generator.qmax, context + ".params");
        generator.mva = requireReal(params, "mva", context + ".params");
        generator.c0  = readRealOrDefault(params, "c0", RealT{0}, context + ".params");
        generator.c1  = readRealOrDefault(params, "c1", RealT{0}, context + ".params");
        generator.c2  = readRealOrDefault(params, "c2", RealT{0}, context + ".params");
        return generator;
      }

      LoadDataT parseLoad(const Json&        buses,
                          const Json&        params,
                          const std::string& id,
                          const std::string& context)
      {
        rejectUnknownFields(buses, {"bus"}, context + ".buses");
        rejectUnknownFields(params,
                            {"pmin", "pmax", "qmin", "qmax"},
                            context + ".params");

        LoadDataT load;
        load.id  = id;
        load.bus = requireUnsignedInteger<IdxT>(buses, "bus", context + ".buses");
        readOptionalReal(params, "pmin", load.pmin, context + ".params");
        readOptionalReal(params, "pmax", load.pmax, context + ".params");
        readOptionalReal(params, "qmin", load.qmin, context + ".params");
        readOptionalReal(params, "qmax", load.qmax, context + ".params");
        return load;
      }

      ShuntDataT parseShunt(const Json&        buses,
                            const Json&        params,
                            const std::string& id,
                            const std::string& context)
      {
        rejectUnknownFields(buses, {"bus"}, context + ".buses");
        rejectUnknownFields(params, {"G", "B"}, context + ".params");

        ShuntDataT shunt;
        shunt.id  = id;
        shunt.bus = requireUnsignedInteger<IdxT>(buses, "bus", context + ".buses");
        shunt.G   = requireReal(params, "G", context + ".params");
        shunt.B   = requireReal(params, "B", context + ".params");
        return shunt;
      }

      DeviceDataT parseDevice(const Json& value, std::size_t position)
      {
        const auto base_context = "devices[" + std::to_string(position) + "]";
        rejectUnknownFields(value, {"class", "buses", "id", "params"}, base_context);

        const auto  device_class = requireString(value, "class", base_context);
        const auto  id           = requireString(value, "id", base_context);
        const auto  context      = base_context + " (" + device_class + " \"" + id + "\")";
        const auto& buses        = requireObjectMember(value, "buses", context);
        const auto& params       = requireObjectMember(value, "params", context);

        if (device_class == "Branch")
        {
          return parseBranch(buses, params, id, context);
        }
        if (device_class == "Generator")
        {
          return parseGenerator(buses, params, id, context);
        }
        if (device_class == "Load")
        {
          return parseLoad(buses, params, id, context);
        }
        if (device_class == "Shunt")
        {
          return parseShunt(buses, params, id, context);
        }

        throwParseError(base_context + ".class",
                        "unsupported device class \"" + device_class + "\"");
      }

      SystemDataT parseSystemDataDocument(const Json& document)
      {
        const std::string context = "root";
        rejectUnknownFields(document, {"header", "params", "buses", "devices"}, context);

        SystemDataT system_data;
        system_data.header = parseHeader(requireObjectMember(document, "header", context));
        system_data.params = parseSystemParameters(
            requireObjectMember(document, "params", context));

        const auto& buses = requireArrayMember(document, "buses", context);
        system_data.buses.reserve(buses.size());
        for (std::size_t i = 0; i < buses.size(); ++i)
        {
          system_data.buses.emplace_back(parseBus(buses[i], i));
        }

        const auto& devices = requireArrayMember(document, "devices", context);
        system_data.devices.reserve(devices.size());
        for (std::size_t i = 0; i < devices.size(); ++i)
        {
          system_data.devices.emplace_back(parseDevice(devices[i], i));
        }
        return system_data;
      }
    } // namespace

    SystemData<> parseSystemData(std::istream& stream)
    {
      try
      {
        return parseSystemDataDocument(Json::parse(stream));
      }
      catch (const Json::exception& error)
      {
        throw std::invalid_argument("OPF system data JSON: " + std::string(error.what()));
      }
    }

    SystemData<> parseSystemData(const std::filesystem::path& file_path)
    {
      std::ifstream stream(file_path);
      if (!stream)
      {
        std::stringstream message;
        message << "Could not open OPF system data file: " << file_path;
        throw std::runtime_error(message.str());
      }

      try
      {
        return parseSystemData(stream);
      }
      catch (const std::invalid_argument& error)
      {
        std::stringstream message;
        message << "Could not parse OPF system data file " << file_path << ": "
                << error.what();
        throw std::invalid_argument(message.str());
      }
    }

  } // namespace OPF
} // namespace GridKit
