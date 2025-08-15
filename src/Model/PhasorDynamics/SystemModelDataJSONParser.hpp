#pragma once

#include <stdexcept>

#include <Model/PhasorDynamics/Bus/BusDataJSONParser.hpp>
#include <Model/PhasorDynamics/ComponentDataJSONParser.hpp>
#include <Model/PhasorDynamics/SignalNode/SignalNodeDataJSONParser.hpp>
#include <Model/PhasorDynamics/SystemModelData.hpp>
#include <nlohmann/json.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using json = nlohmann::json;

    /// JSON parser function implementation for the `SystemModelData` type
    ///
    /// See the `README.md` in `src/Model/PhasorDynamics` for more information
    template <typename RealT = double, typename IdxT = size_t>
    void from_json(const json& j, SystemModelData<RealT, IdxT>& sm)
    {
      auto header = j.at("header");

      if (header.contains("format_version"))
      {
        header.at("format_version").get_to(sm.format_version);
      }

      if (header.contains("format_revision"))
      {
        header.at("format_revision").get_to(sm.format_revision);
      }

      header.at("case_name").get_to(sm.case_name);

      if (header.contains("case_date_time"))
      {
        header.at("case_date_time").get_to(sm.case_date_time);
      }

      header.at("case_description").get_to(sm.case_description);
      header.at("case_comments").get_to(sm.case_comments);
      header.at("freq_base").get_to(sm.freq_base);
      header.at("va_base").get_to(sm.va_base);

      /// @todo Give signal nodes their own array!!!
      /// Modify JSON format accordingly

      /// Gets all electrical buses
      j.at("buses").get_to(sm.bus);

      /// Gets all signal nodes (allows for systems without signals)
      if (j.contains("signals"))
      {
        j.at("signals").get_to(sm.signal);
      }

      /// Gets all components
      /// @todo So far handling only branches, Genrous, and bus faults
      /// For 2-bus Tgov1, for example, we need to add support for Tgov1
      /// governor and support for signal nodes.
      for (auto& raw_component : j.at("devices"))
      {
        auto kind = raw_component.at("class");
        if (kind == "branch")
        {
          typename SystemModelData<RealT, IdxT>::BranchDataT branch;
          raw_component.get_to(branch);
          sm.branch.push_back(branch);
        }
        else if (kind == "GENROU")
        {
          typename SystemModelData<RealT, IdxT>::GenrouDataT genrou;
          raw_component.get_to(genrou);
          sm.genrou.push_back(genrou);
        }
        else if (kind == "bus_fault")
        {
          typename SystemModelData<RealT, IdxT>::BusFaultDataT bus_fault;
          raw_component.get_to(bus_fault);
          sm.bus_fault.push_back(bus_fault);
        }
        else
        {
          throw std::runtime_error("Invalid device class");
        }
      }
    }
  } // namespace PhasorDynamics
} // namespace GridKit
