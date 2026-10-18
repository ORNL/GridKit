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
    /**
     * @brief All data needed to build an optimal power flow system model
     */
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

      std::vector<BusDataT>       bus;       ///< Buses within the model
      std::vector<BranchDataT>    branch;    ///< Branches within the model
      std::vector<GeneratorDataT> generator; ///< Generators within the model
      std::vector<LoadDataT>      load;      ///< Loads within the model
    };

    /**
     * @brief Generate system model data from a JSON input file
     *
     * @throw std::invalid_argument for malformed or unsupported input,
     * naming the bus or device
     */
    SystemModelData<double, size_t> parseSystemModelData(std::istream& stream);
    SystemModelData<double, size_t> parseSystemModelData(const std::filesystem::path& file);
  } // namespace OptimalPowerFlow
} // namespace GridKit
