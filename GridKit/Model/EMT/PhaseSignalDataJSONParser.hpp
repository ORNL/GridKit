#pragma once

#include <GridKit/Model/EMT/ComponentDataJSONParser.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Expand a documented vector port into the scalar signal map.
    template <size_t N = 3>
    inline void expandPhasePort(json& j, const char* direction, const char* port, const std::array<const char*, N>& phases)
    {
      if (!j.contains(direction) || !j.at(direction).contains(port))
      {
        return;
      }
      auto&      ports  = j.at(direction);
      const auto values = ports.at(port).get<std::vector<std::string>>();
      if (values.size() != N)
      {
        throw std::invalid_argument(std::string(port) + " requires " + std::to_string(N) + " signal IDs");
      }
      for (size_t n = 0; n < N; ++n)
      {
        if (values[n].empty() || ports.contains(phases[n]))
        {
          throw std::invalid_argument(std::string(port) + " has an empty or duplicate phase mapping");
        }
        ports[phases[n]] = values[n];
      }
      ports.erase(port);
    }

    /// Vector monitors use one scalar column per phase in every output format.
    template <size_t N = 3>
    inline void expandPhaseMonitor(json& j, const char* name, const std::array<const char*, N>& phases)
    {
      if (!j.contains("mon"))
      {
        return;
      }
      if (!j.at("mon").is_array())
        throw std::invalid_argument("mon must be an array");
      auto monitors = json::array();
      for (const auto& monitor : j.at("mon"))
      {
        if (monitor == name)
        {
          for (const auto* phase : phases)
          {
            monitors.push_back(phase);
          }
        }
        else
        {
          monitors.push_back(monitor);
        }
      }
      j["mon"] = std::move(monitors);
    }
  } // namespace EMT
} // namespace GridKit
