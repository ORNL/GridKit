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
#include <GridKit/Model/OptimalPowerFlow/Controller/REECB/ReecbData.hpp>
#include <GridKit/Model/OptimalPowerFlow/Controller/REPCA/RepcaData.hpp>
#include <GridKit/Model/OptimalPowerFlow/Converter/REGCA/RegcaData.hpp>
#include <GridKit/Model/OptimalPowerFlow/Exciter/ESDC1A/Esdc1aData.hpp>
#include <GridKit/Model/OptimalPowerFlow/Exciter/IEEET1/Ieeet1Data.hpp>
#include <GridKit/Model/OptimalPowerFlow/Exciter/SEXS-PTI/SexsPtiData.hpp>
#include <GridKit/Model/OptimalPowerFlow/Generator/GeneratorData.hpp>
#include <GridKit/Model/OptimalPowerFlow/Governor/GASTPTI/GastPtiData.hpp>
#include <GridKit/Model/OptimalPowerFlow/Governor/HYGOV/HygovData.hpp>
#include <GridKit/Model/OptimalPowerFlow/Governor/Tgov1/Tgov1Data.hpp>
#include <GridKit/Model/OptimalPowerFlow/Load/LoadData.hpp>
#include <GridKit/Model/OptimalPowerFlow/SynchronousMachine/GENROU/GenrouData.hpp>
#include <GridKit/Model/OptimalPowerFlow/SynchronousMachine/GENSAL/GensalData.hpp>
#include <GridKit/Model/OptimalPowerFlow/SynchronousMachine/GenClassical/GenClassicalData.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /// Data of an optimal power flow system
    template <typename real_type = double, typename index_type = size_t>
    struct SystemModelData
    {
      using RealT             = real_type;
      using IdxT              = index_type;
      using BusDataT          = BusData<RealT, IdxT>;
      using BranchDataT       = BranchData<RealT, IdxT>;
      using GeneratorDataT    = GeneratorData<RealT, IdxT>;
      using LoadDataT         = LoadData<RealT, IdxT>;
      using GenrouDataT       = GenrouData<RealT, IdxT>;
      using GensalDataT       = GensalData<RealT, IdxT>;
      using GenClassicalDataT = GenClassicalData<RealT, IdxT>;
      using Esdc1aDataT       = Exciter::Esdc1aData<RealT, IdxT>;
      using Ieeet1DataT       = Exciter::Ieeet1Data<RealT, IdxT>;
      using SexsPtiDataT      = Exciter::SexsPtiData<RealT, IdxT>;
      using GastPtiDataT      = Governor::GastPtiData<RealT, IdxT>;
      using HygovDataT        = Governor::HygovData<RealT, IdxT>;
      using Tgov1DataT        = Governor::Tgov1Data<RealT, IdxT>;
      using RegcaDataT        = Converter::RegcaData<RealT, IdxT>;
      using ReecbDataT        = Controller::ReecbData<RealT, IdxT>;
      using RepcaDataT        = Controller::RepcaData<RealT, IdxT>;

      RealT va_base{100.0e6}; ///< System power base in VA

      std::vector<BusDataT>          bus;
      std::vector<IdxT>              signal; ///< Signal ids
      std::vector<BranchDataT>       branch;
      std::vector<GeneratorDataT>    generator;
      std::vector<LoadDataT>         load;
      std::vector<GenrouDataT>       genrou;
      std::vector<GensalDataT>       gensal;
      std::vector<GenClassicalDataT> gen_classical;
      std::vector<Esdc1aDataT>       esdc1a;
      std::vector<Ieeet1DataT>       ieeet1;
      std::vector<SexsPtiDataT>      sexs_pti;
      std::vector<GastPtiDataT>      gast_pti;
      std::vector<HygovDataT>        hygov;
      std::vector<Tgov1DataT>        tgov1;
      std::vector<RegcaDataT>        regca;
      std::vector<ReecbDataT>        reecb;
      std::vector<RepcaDataT>        repca;
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
