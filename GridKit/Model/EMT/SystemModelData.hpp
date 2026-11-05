/**
 * @file SystemModelData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Data of an EMT system.
 */

#pragma once

#include <cstddef>
#include <filesystem>
#include <istream>
#include <map>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include <GridKit/Constants.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Bus or device of a case
    template <typename index_type>
    struct DeviceData
    {
      std::string                       device_class;
      std::string                       id;
      nlohmann::json                    params;
      std::map<std::string, index_type> ports; ///< Bus number at each port
      std::vector<std::string>          mon;   ///< Monitored internal variables
    };

    template <typename real_type = double, typename index_type = size_t>
    struct SystemModelData
    {
      using RealT           = real_type;
      using IdxT            = index_type;
      using MonitorSinkSpec = Model::VariableMonitorBase::SinkSpec;

      RealT va_base{ZERO<RealT>}; ///< System power base [VA]

      std::vector<DeviceData<IdxT>> buses; ///< Each with its number as id
      std::vector<DeviceData<IdxT>> devices;

      std::vector<MonitorSinkSpec> monitor_sink;
    };

    /// @throw nlohmann::json::exception for a missing or mistyped value
    SystemModelData<double, size_t> parseSystemModelData(std::istream& stream);
    SystemModelData<double, size_t> parseSystemModelData(const std::filesystem::path& file);
  } // namespace EMT
} // namespace GridKit
