#pragma once

#include <optional>

#include <Model/PhasorDynamics/Branch/BranchData.hpp>
#include <Model/PhasorDynamics/Bus/BusData.hpp>
#include <Model/PhasorDynamics/BusFault/BusFaultData.hpp>
#include <Model/PhasorDynamics/Governor/Tgov1/Tgov1Data.hpp>
#include <Model/PhasorDynamics/Load/LoadData.hpp>
#include <Model/PhasorDynamics/SynchronousMachine/GENROUwS/GenrouData.hpp>
#include <nlohmann/json.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using json = nlohmann::json;

    /// A structure containing all of the data needed to reproduce a
    /// system model
    ///
    /// In particular, this structure is modeled after the grid dynamics
    /// case format, which is described within `INPUT_FORMAT.md`
    template <typename RealT = double, typename IdxT = size_t>
    struct SystemModelData
    {
      using BranchDataT   = BranchData<RealT, IdxT>;
      using BusDataT      = BusData<RealT, IdxT>;
      using BusFaultDataT = BusFaultData<RealT, IdxT>;
      using Tgov1DataT    = Governor::Tgov1Data<RealT, IdxT>;
      using GenrouDataT   = GenrouData<RealT, IdxT>;
      using LoadDataT     = LoadData<RealT, IdxT>;

      /// The version of the grid dynamics case format this system model was
      /// parsed from
      ///
      /// If not relevant (i.e. if not working in the context of the case
      /// format), this will not contain a value
      std::optional<unsigned short> format_version;

      /// The revision of the grid dynamics case format this system model was
      /// parsed from
      ///
      /// If not relevant (i.e. not working in the context of the case
      /// format), this will not contain a value
      std::optional<unsigned short> format_revision;

      /// A name for the case this model describes
      std::string case_name;

      /// A date and time relevant to the case being described by this model
      ///
      /// TODO: this should probably use the chrono header instead of just
      ///       being a string
      std::optional<std::string> case_date_time;

      /// A description of the case this model describes
      std::string case_description;

      /// Additional comments about the case being described by this model
      std::string case_comments;

      RealT freq_base; ///< System frequency base in Hz
      RealT va_base;   ///< System power base in VA

      std::vector<BusDataT>      bus;       ///< Buses within the model
      std::vector<BranchDataT>   branch;    ///< Branches within the model
      std::vector<BusFaultDataT> bus_fault; ///< Bus faults within the model
      std::vector<GenrouDataT>   genrou;    ///< GENROU instances within the model
      std::vector<LoadDataT>     load;      ///< Loads within the model
      std::vector<Tgov1DataT>    gov;       ///< Governors within the model
    };

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
      j.at("buses").get_to(sm.bus);

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
          throw "Invalid device class";
        }
      }
    }
  } // namespace PhasorDynamics
} // namespace GridKit
