#pragma once

#include <Model/PhasorDynamics/Branch/BranchData.hpp>
#include <Model/PhasorDynamics/Bus/BusData.hpp>
#include <Model/PhasorDynamics/BusFault/BusFaultData.hpp>
#include <Model/PhasorDynamics/Load/LoadData.hpp>
#include <Model/PhasorDynamics/SynchronousMachine/GENROUwS/GenrouData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename RealT = double, typename IdxT = size_t>
    struct SystemModelData
    {
      using BranchDataT   = BranchData<RealT, IdxT>;
      using BusDataT      = BusData<RealT, IdxT>;
      using BusFaultDataT = BusFaultData<RealT, IdxT>;
      using GenrouDataT   = GenrouData<RealT, IdxT>;
      using LoadDataT     = LoadData<RealT, IdxT>;

      // Header info
      unsigned short format_version;
      unsigned short format_revision;
      std::string    case_name;
      std::string    case_date_time;
      std::string    case_description;
      std::string    case_comments;
      RealT          freq_base;
      RealT          va_base;

      // Bus data
      std::vector<BusDataT> bus;

      // Component data
      std::vector<BranchDataT>   branch;
      std::vector<BusFaultDataT> bus_fault;
      std::vector<GenrouDataT>   genrou;
      std::vector<LoadDataT>     load;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
