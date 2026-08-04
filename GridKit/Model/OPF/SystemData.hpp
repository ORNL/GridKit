#pragma once

#include <cstddef>
#include <filesystem>
#include <istream>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include <GridKit/Model/OPF/Branch/BranchData.hpp>
#include <GridKit/Model/OPF/Bus/BusData.hpp>
#include <GridKit/Model/OPF/Generator/GeneratorData.hpp>
#include <GridKit/Model/OPF/Load/LoadData.hpp>
#include <GridKit/Model/OPF/Shunt/ShuntData.hpp>

namespace GridKit
{
  namespace OPF
  {
    struct CaseHeader
    {
      unsigned int               format_version{};
      unsigned int               format_revision{};
      std::string                case_name;
      std::optional<std::string> case_date_time;
      std::optional<std::string> case_description;
      std::optional<std::string> case_comments;
    };

    template <typename real_type>
    struct SystemParameters
    {
      using RealT = real_type;

      RealT freq_base{};
      RealT va_base{};
    };

    template <typename real_type = double, typename index_type = std::size_t>
    struct SystemData
    {
      using RealT          = real_type;
      using IdxT           = index_type;
      using BusDataT       = BusData<RealT, IdxT>;
      using BranchDataT    = BranchData<RealT, IdxT>;
      using GeneratorDataT = GeneratorData<RealT, IdxT>;
      using LoadDataT      = LoadData<RealT, IdxT>;
      using ShuntDataT     = ShuntData<RealT, IdxT>;
      using DeviceDataT    = std::variant<BranchDataT,
                                          GeneratorDataT,
                                          LoadDataT,
                                          ShuntDataT>;

      CaseHeader               header;
      SystemParameters<RealT>  params;
      std::vector<BusDataT>    buses;
      std::vector<DeviceDataT> devices;
    };

    SystemData<> parseSystemData(std::istream& stream);

    SystemData<> parseSystemData(const std::filesystem::path& file_path);

  } // namespace OPF
} // namespace GridKit
