#pragma once

#include <filesystem>
#include <istream>
#include <optional>
#include <string>
#include <vector>

#include <GridKit/Model/PhasorDynamics/Branch/BranchData.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusData.hpp>
#include <GridKit/Model/PhasorDynamics/BusFault/BusFaultData.hpp>
#include <GridKit/Model/PhasorDynamics/BusToSignalAdapter/BusToSignalAdapterData.hpp>
#include <GridKit/Model/PhasorDynamics/Controller/REPCA/RepcaData.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REGCA/RegcaData.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/ESDC1A/Esdc1aData.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/IEEET1/Ieeet1Data.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/SEXS-PTI/SexsPtiData.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/GASTPTI/GastPtiData.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/HYGOV/HygovData.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/Tgov1/Tgov1Data.hpp>
#include <GridKit/Model/PhasorDynamics/Load/LoadZ/LoadZData.hpp>
#include <GridKit/Model/PhasorDynamics/Load/LoadZIP/LoadZIPData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalSource/ConstantSignalSourceData.hpp>
#include <GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/IeeestData.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENROU/GenrouData.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENSAL/GensalData.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GenClassical/GenClassicalData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// A structure containing all of the data needed to reproduce a
    /// system model
    ///
    /// In particular, this structure is modeled after the grid dynamics
    /// case format, which is described within `INPUT_FORMAT.md`
    template <typename real_type = double, typename index_type = size_t>
    struct SystemModelData
    {
      using RealT                   = real_type;
      using IdxT                    = index_type;
      using BranchDataT             = BranchData<RealT, IdxT>;
      using BusDataT                = BusData<RealT, IdxT>;
      using BusToSignalAdapterDataT = BusToSignalAdapterData<RealT, IdxT>;
      using BusFaultDataT           = BusFaultData<RealT, IdxT>;
      using RegcaDataT              = Converter::RegcaData<RealT, IdxT>;
      using RepcaDataT              = Controller::RepcaData<RealT, IdxT>;
      using Tgov1DataT              = Governor::Tgov1Data<RealT, IdxT>;
      using Esdc1aDataT             = Exciter::Esdc1aData<RealT, IdxT>;
      using GastPtiDataT            = Governor::GastPtiData<RealT, IdxT>;
      using HygovDataT              = Governor::HygovData<RealT, IdxT>;
      using Ieeet1DataT             = Exciter::Ieeet1Data<RealT, IdxT>;
      using SexsPtiDataT            = Exciter::SexsPtiData<RealT, IdxT>;
      using IeeestDataT             = Stabilizer::IeeestData<RealT, IdxT>;
      using GenrouDataT             = GenrouData<RealT, IdxT>;
      using GensalDataT             = GensalData<RealT, IdxT>;
      using GenClassicalDataT       = GenClassicalData<RealT, IdxT>;
      using LoadZDataT              = LoadZData<RealT, IdxT>;
      using LoadZIPDataT            = LoadZIPData<RealT, IdxT>;
      using ConstantSourceT         = ConstantSignalSourceData<RealT, IdxT>;
      using SignalDataT             = SignalNodeData<RealT, IdxT>;
      using MonitorSinkSpec         = Model::VariableMonitorBase::SinkSpec;

      /// The version of the grid dynamics case format this system model was
      /// parsed from
      ///
      /// If not relevant (i.e. if not working in the context of the case
      /// format), this will not contain a value
      std::optional<double> format_version;

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

      RealT freq_base{60.0};  ///< System frequency base in Hz
      RealT va_base{100.0e6}; ///< System power base in VA

      /// @todo Create an enum identifying all available component models and
      /// consolidate components in a single container.
      /// - Convert string to enum
      /// - Associate component type to its corresponding enum
      /// - Consolidate components to allow writing to them using the enum as the argument
      std::vector<BusDataT>                bus;             ///< Buses within the model
      std::vector<BusToSignalAdapterDataT> adapter;         ///< bus-to-signal adapters within the model
      std::vector<BranchDataT>             branch;          ///< Branches within the model
      std::vector<BusFaultDataT>           bus_fault;       ///< Bus faults within the model
      std::vector<RegcaDataT>              regca;           ///< REGCA converter instances within the model
      std::vector<RepcaDataT>              repca;           ///< REPCA plant controllers within the model
      std::vector<GenrouDataT>             genrou;          ///< GENROU instances within the model
      std::vector<GensalDataT>             gensal;          ///< GENSAL instances within the model
      std::vector<GenClassicalDataT>       genclassical;    ///< Classical generator instances within the model
      std::vector<LoadZDataT>              loadz;           ///< LoadZ instances within the model
      std::vector<LoadZIPDataT>            loadzip;         ///< LoadZIP instances within the model
      std::vector<Tgov1DataT>              gov;             ///< Governors within the model
      std::vector<Esdc1aDataT>             esdc1a;          ///< ESDC1A exciters within the model
      std::vector<GastPtiDataT>            gastpti;         ///< GASTPTI governors within the model
      std::vector<HygovDataT>              hygov;           ///< HYGOV governors within the model
      std::vector<Ieeet1DataT>             exciter;         ///< Exciters within the model
      std::vector<SexsPtiDataT>            sexspti;         ///< SEXS-PTI exciters within the model
      std::vector<IeeestDataT>             stabilizer;      ///< Stabilizers within the model
      std::vector<ConstantSourceT>         constant_source; ///< Constant signal sources within the model
      std::vector<SignalDataT>             signal;          ///< Signal nodes

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
  } // namespace PhasorDynamics
} // namespace GridKit
