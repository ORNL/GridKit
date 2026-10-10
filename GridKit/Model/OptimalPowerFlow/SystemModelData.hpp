/**
 * @file SystemModelData.hpp
 * @brief Data of an optimal power flow system model.
 */

#pragma once

#include <cstddef>
#include <vector>

#include <GridKit/Model/OptimalPowerFlow/Branch/BranchData.hpp>
#include <GridKit/Model/OptimalPowerFlow/Bus/BusData.hpp>
#include <GridKit/Model/OptimalPowerFlow/Generator/GeneratorData.hpp>
#include <GridKit/Model/OptimalPowerFlow/Load/LoadData.hpp>
#include <GridKit/Model/OptimalPowerFlow/MatpowerData.hpp>
#include <GridKit/Model/OptimalPowerFlow/Shunt/ShuntData.hpp>

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
      using ShuntDataT     = ShuntData<RealT, IdxT>;

      RealT va_base{100.0e6}; ///< System power base in VA

      std::vector<BusDataT>       bus;       ///< Buses within the model
      std::vector<BranchDataT>    branch;    ///< Branches within the model
      std::vector<GeneratorDataT> generator; ///< Generators within the model
      std::vector<LoadDataT>      load;      ///< Loads within the model
      std::vector<ShuntDataT>     shunt;     ///< Shunts within the model
    };

    /**
     * @brief Set bus voltage limits, branch ratings, and generator limits and
     * costs from a MATPOWER case of the same network
     *
     * Buses match by number. In-service MATPOWER branches match the branches
     * between the same buses in order, and in-service MATPOWER generators the
     * generators at the same bus. MATPOWER generators at a bus without
     * generators are static injections that the loads carry, and are skipped.
     * Values in MW, Mvar, and MVA become system base values.
     */
    void applyMatpowerData(SystemModelData<double, size_t>& data, const MatpowerData& matpower);
  } // namespace OptimalPowerFlow
} // namespace GridKit
