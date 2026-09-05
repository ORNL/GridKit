#pragma once

#include <GridKit/Model/EMT/ComponentDataJSONParser.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Expand a documented three-phase port into the scalar signal map.
    inline void expandPhasePort(json& j, const char* direction, const char* port, const std::array<const char*, 3>& phases)
    {
      if (!j.contains(direction) || !j.at(direction).contains(port))
      {
        return;
      }
      auto&      ports  = j.at(direction);
      const auto values = ports.at(port).get<std::vector<std::string>>();
      if (values.size() != 3)
      {
        throw std::invalid_argument(std::string(port) + " requires three signal IDs");
      }
      for (size_t n = 0; n < 3; ++n)
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
    inline void expandPhaseMonitor(json& j, const char* name, const std::array<const char*, 3>& phases)
    {
      if (!j.contains("mon"))
      {
        return;
      }
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
