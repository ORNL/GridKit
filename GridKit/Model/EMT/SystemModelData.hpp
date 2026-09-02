#pragma once

#include <filesystem>
#include <istream>
#include <optional>
#include <string>
#include <vector>

#include <GridKit/Model/EMT/Bus/BusData.hpp>
#include <GridKit/Model/EMT/Component/Line/LineLumped/LineLumpedData.hpp>
#include <GridKit/Model/EMT/Component/Load/LoadZ/LoadZData.hpp>
#include <GridKit/Model/EMT/Component/Machine/SynchronousMachine/SynchronousMachineData.hpp>
#include <GridKit/Model/EMT/Component/Source/VoltageSource/VoltageSourceData.hpp>
#include <GridKit/Model/EMT/Component/Switch/SwitchData.hpp>
#include <GridKit/Model/EMT/Governor/Tgov1/Tgov1Data.hpp>
#include <GridKit/Model/EMT/SignalNode/SignalNodeData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// A structure containing all of the data needed to reproduce a
    /// system model
    ///
    /// In particular, this structure is modeled after the EMT case format,
    /// which is described within `INPUT_FORMAT.md`
    template <typename real_type = double, typename index_type = size_t>
    struct SystemModelData
    {
      using RealT                   = real_type;
      using IdxT                    = index_type;
      using BusDataT                = BusData<RealT, IdxT>;
      using LineLumpedDataT         = LineLumpedData<RealT, IdxT>;
      using LoadZDataT              = LoadZData<RealT, IdxT>;
      using SynchronousMachineDataT = SynchronousMachineData<RealT, IdxT>;
      using VoltageSourceDataT      = VoltageSourceData<RealT, IdxT>;
      using SwitchDataT             = SwitchData<RealT, IdxT>;
      using Tgov1DataT              = Governor::Tgov1Data<RealT, IdxT>;
      using SignalDataT             = SignalNodeData<RealT, IdxT>;
      using MonitorSinkSpec         = Model::VariableMonitorBase::SinkSpec;

      /// The version of the EMT case format this system model was parsed from
      ///
      /// If not relevant (i.e. if not working in the context of the case
      /// format), this will not contain a value
      std::optional<double> format_version;

      /// The revision of the EMT case format this system model was parsed from
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

      std::vector<BusDataT>                bus;                 ///< Buses within the model
      std::vector<VoltageSourceDataT>      voltage_source;      ///< Voltage sources within the model
      std::vector<SynchronousMachineDataT> synchronous_machine; ///< Synchronous machines within the model
      std::vector<LineLumpedDataT>         line_lumped;         ///< Lumped lines within the model
      std::vector<LoadZDataT>              loadz;               ///< LoadZ instances within the model
      std::vector<SwitchDataT>             sw;                  ///< Switches within the model
      std::vector<Tgov1DataT>              gov;                 ///< Governors within the model
      std::vector<SignalDataT>             signal;              ///< Signal nodes

      /// Monitor sink specs
      std::vector<MonitorSinkSpec> monitor_sink;
    };

    /**
     * @brief Generate system model data from a JSON input file
     */
    SystemModelData<double, size_t> parseSystemModelData(std::istream& stream);
    SystemModelData<double, size_t> parseSystemModelData(std::istream&& stream);
    SystemModelData<double, size_t> parseSystemModelData(const std::filesystem::path& filePath);
    SystemModelData<double, size_t> parseSystemModelData(const std::string& fileName);
  } // namespace EMT
} // namespace GridKit
