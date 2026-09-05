#pragma once

#include <filesystem>
#include <istream>
#include <optional>
#include <string>
#include <vector>

#include <GridKit/Model/EMT/ContainerData.hpp>
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
    struct SystemModelData : public ContainerData<real_type, index_type>
    {
      using RealT = real_type;
      using IdxT  = index_type;

      using MonitorSinkSpec = Model::VariableMonitorBase::SinkSpec;

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
