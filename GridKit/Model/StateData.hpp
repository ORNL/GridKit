#pragma once

#include <filesystem>
#include <istream>
#include <map>
#include <optional>
#include <ostream>
#include <string>

namespace GridKit
{
  namespace Model
  {
    /// Metadata associated with a model state artifact.
    struct HeaderState
    {
      std::optional<unsigned int> version;     ///< Descriptive value from the version field
      std::optional<double>       time;        ///< Model measurement time
      std::optional<std::string>  created;     ///< Wall-clock creation time
      std::optional<std::string>  description; ///< Description of the state
    };

    /// Current injected by a device into a bus.
    struct InjectionState
    {
      std::optional<double> ir; ///< Real current
      std::optional<double> ii; ///< Imaginary current
      std::optional<double> ia; ///< Instantaneous phase a current
      std::optional<double> ib; ///< Instantaneous phase b current
      std::optional<double> ic; ///< Instantaneous phase c current
    };

    /// State associated with a bus.
    struct BusState
    {
      std::optional<double>                 vr; ///< Real voltage
      std::optional<double>                 vi; ///< Imaginary voltage
      std::optional<double>                 va; ///< Instantaneous phase a voltage
      std::optional<double>                 vb; ///< Instantaneous phase b voltage
      std::optional<double>                 vc; ///< Instantaneous phase c voltage
      std::map<std::string, InjectionState> injections;
    };

    /// State associated with a device.
    struct DeviceState
    {
      std::optional<bool>   active;
      std::optional<bool>   online;
      std::optional<bool>   open;
      std::optional<double> p;
      std::optional<double> q;
      std::optional<double> tap;
      std::optional<double> phase;
    };

    /// Portable model state parsed from or written to a state artifact.
    struct StateData
    {
      std::optional<HeaderState>         header;
      std::map<std::string, BusState>    buses;
      std::map<std::string, DeviceState> devices;
    };

    /// Parse model state from a JSON stream.
    StateData parseStateData(std::istream& stream);

    /// Parse model state from a JSON file.
    StateData parseStateData(const std::filesystem::path& file_path);

    /// Write model state as JSON to a stream.
    void writeStateData(const StateData& state, std::ostream& stream);

    /// Write model state as JSON to a file.
    void writeStateData(const StateData&             state,
                        const std::filesystem::path& file_path);
  } // namespace Model
} // namespace GridKit
