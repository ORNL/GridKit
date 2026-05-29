#pragma once

#include <sstream>
#include <stdexcept>

#include <nlohmann/json.hpp>

#include <GridKit/Model/PhasorDynamics/Bus/BusDataJSONParser.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentDataJSONParser.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeDataJSONParser.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using json          = nlohmann::json;
    using Log           = ::GridKit::Utilities::Logger;
    using MonitorFormat = ::GridKit::Model::VariableMonitorFormat;

    /// JSON parser function implementation for the `SystemModelData` type
    ///
    /// See the `README.md` in `GridKit/Model/PhasorDynamics` for more information
    template <typename RealT = double, typename IdxT = size_t>
    void from_json(const json& j, SystemModelData<RealT, IdxT>& sm)
    {
      auto enum_parse = []<typename EnumT, typename KeyT>(EnumT, KeyT&& key)
      {
        return magic_enum::enum_cast<EnumT>(key, magic_enum::case_insensitive);
      };

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

      if (j.contains("monitors"))
      {
        for (auto&& raw_mon : j.at("monitors"))
        {
          auto file_name = raw_mon.value("file_name", std::string{});
          auto fmt_str   = raw_mon.at("format").get<std::string>();
          auto format    = enum_parse(MonitorFormat{}, fmt_str);
          auto delim     = raw_mon.value("delim", std::string(","));
          if (format.has_value())
          {
            sm.monitor_sink.emplace_back(format.value(), file_name, delim);
          }
          else
          {
            Log::error() << "\n\tInvalid monitor format: \"" << fmt_str << "\"."
                         << "\n\tSee the \"monitors\" list in your JSON file."
                         << std::endl;
          }
        }
      }

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
        auto kind = raw_component.at("class").get<std::string>();
        if (kind == "BusToSignalAdapter")
        {
          typename SystemModelData<RealT, IdxT>::BusToSignalAdapterDataT adapter;
          raw_component.get_to(adapter);
          sm.adapter.push_back(adapter);
        }
        if (kind == "Branch")
        {
          typename SystemModelData<RealT, IdxT>::BranchDataT branch;
          raw_component.get_to(branch);
          sm.branch.push_back(branch);
        }
        else if (kind == "Genrou")
        {
          typename SystemModelData<RealT, IdxT>::GenrouDataT genrou;
          raw_component.get_to(genrou);
          sm.genrou.push_back(genrou);
        }
        else if (kind == "Gensal")
        {
          typename SystemModelData<RealT, IdxT>::GensalDataT gensal;
          raw_component.get_to(gensal);
          sm.gensal.push_back(gensal);
        }
        else if (kind == "GenClassical")
        {
          typename SystemModelData<RealT, IdxT>::GenClassicalDataT gen_classical;
          raw_component.get_to(gen_classical);
          sm.genclassical.push_back(gen_classical);
        }
        else if (kind == "Load")
        {
          typename SystemModelData<RealT, IdxT>::LoadDataT load;
          raw_component.get_to(load);
          sm.load.push_back(load);
        }
        else if (kind == "LoadZIP")
        {
          typename SystemModelData<RealT, IdxT>::LoadZIPDataT loadzip;
          raw_component.get_to(loadzip);
          sm.loadzip.push_back(loadzip);
        }
        else if (kind == "Tgov1")
        {
          typename SystemModelData<RealT, IdxT>::Tgov1DataT gov;
          raw_component.get_to(gov);
          sm.gov.push_back(gov);
        }
        else if (kind == "Hygov")
        {
          typename SystemModelData<RealT, IdxT>::HygovDataT gov;
          raw_component.get_to(gov);
          sm.hygov.push_back(gov);
        }
        else if (kind == "Ieeet1")
        {
          typename SystemModelData<RealT, IdxT>::Ieeet1DataT exciter;
          raw_component.get_to(exciter);
          sm.exciter.push_back(exciter);
        }
        else if (kind == "SexsPti")
        {
          typename SystemModelData<RealT, IdxT>::SexsPtiDataT exciter;
          raw_component.get_to(exciter);
          sm.sexspti.push_back(exciter);
        }
        else if (kind == "Ieeest")
        {
          typename SystemModelData<RealT, IdxT>::IeeestDataT stabilizer;
          raw_component.get_to(stabilizer);
          sm.stabilizer.push_back(stabilizer);
        }
        else if (kind == "BusFault")
        {
          typename SystemModelData<RealT, IdxT>::BusFaultDataT bus_fault;
          raw_component.get_to(bus_fault);
          sm.bus_fault.push_back(bus_fault);
        }
        else
        {
          Log::error() << "\n\tInvalid device class: \"" << kind << "\". "
                       << "\n\tSee the \"devices\" list in your JSON file."
                       << std::endl;
        }
      }
    }
  } // namespace PhasorDynamics
} // namespace GridKit
