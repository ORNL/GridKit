#pragma once

#include <optional>

#include <Model/PhasorDynamics/Branch/BranchData.hpp>
#include <Model/PhasorDynamics/Bus/BusData.hpp>
#include <Model/PhasorDynamics/BusFault/BusFaultData.hpp>
#include <Model/PhasorDynamics/Load/LoadData.hpp>
#include <Model/PhasorDynamics/SynchronousMachine/GENROUwS/GenrouData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
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

      /// System frequency base in Hz
      RealT freq_base;

      /// System power base in VA
      RealT va_base;

      /// Buses within the model
      std::vector<BusDataT> bus;

      // Branches within the model
      std::vector<BranchDataT> branch;

      /// Bus faults within the model
      std::vector<BusFaultDataT> bus_fault;

      /// GENROU synchronous machine model instances within the model
      std::vector<GenrouDataT> genrou;

      /// Loads within the model
      std::vector<LoadDataT> load;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
