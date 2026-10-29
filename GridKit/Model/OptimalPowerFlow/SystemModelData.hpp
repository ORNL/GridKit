/**
 * @file SystemModelData.hpp
 * @brief Data of an optimal power flow system model.
 */

#pragma once

#include <cstddef>
#include <filesystem>
#include <istream>
#include <vector>

#include <GridKit/Model/OptimalPowerFlow/Branch/BranchData.hpp>
#include <GridKit/Model/OptimalPowerFlow/Bus/BusData.hpp>
#include <GridKit/Model/OptimalPowerFlow/Generator/GeneratorData.hpp>
#include <GridKit/Model/OptimalPowerFlow/Load/LoadData.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /// Data of an optimal power flow system
    template <typename real_type = double, typename index_type = size_t>
    struct SystemModelData
    {
      using RealT          = real_type;
      using IdxT           = index_type;
      using BusDataT       = BusData<RealT, IdxT>;
      using BranchDataT    = BranchData<RealT, IdxT>;
      using GeneratorDataT = GeneratorData<RealT, IdxT>;
      using LoadDataT      = LoadData<RealT, IdxT>;

      RealT va_base{100.0e6}; ///< System power base in VA

      std::vector<BusDataT>       bus;
      std::vector<BranchDataT>    branch;
      std::vector<GeneratorDataT> generator;
      std::vector<LoadDataT>      load;
    };

    /**
     * @brief System data from a JSON case
     *
     * @throw nlohmann::json::exception for a missing or mistyped value, and
     * std::invalid_argument for an unknown device class, parameter, or port
     */
    SystemModelData<double, size_t> parseSystemModelData(std::istream& stream);
    SystemModelData<double, size_t> parseSystemModelData(const std::filesystem::path& file);
  } // namespace OptimalPowerFlow
} // namespace GridKit
