#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <variant>
#include <vector>

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

    enum class BusClass
    {
      BUS,
      SLACK
    };

    template <typename real_type, typename index_type>
    struct BusData
    {
      using RealT = real_type;
      using IdxT  = index_type;

      BusClass             bus_class{BusClass::BUS};
      std::string          id;
      IdxT                 number{};
      RealT                kv{};
      std::optional<RealT> vmin;
      std::optional<RealT> vmax;
    };

    template <typename real_type, typename index_type>
    struct BranchData
    {
      using RealT = real_type;
      using IdxT  = index_type;

      std::string          id;
      IdxT                 from{};
      IdxT                 to{};
      RealT                R{};
      RealT                X{};
      RealT                G{};
      RealT                B{};
      std::optional<RealT> smax;
    };

    template <typename real_type, typename index_type>
    struct GeneratorData
    {
      using RealT = real_type;
      using IdxT  = index_type;

      std::string          id;
      IdxT                 bus{};
      std::optional<RealT> pmin;
      std::optional<RealT> pmax;
      std::optional<RealT> qmin;
      std::optional<RealT> qmax;
      RealT                mva{};
      RealT                c0{};
      RealT                c1{};
      RealT                c2{};
    };

    template <typename real_type, typename index_type>
    struct LoadData
    {
      using RealT = real_type;
      using IdxT  = index_type;

      std::string          id;
      IdxT                 bus{};
      std::optional<RealT> pmin;
      std::optional<RealT> pmax;
      std::optional<RealT> qmin;
      std::optional<RealT> qmax;
    };

    template <typename real_type, typename index_type>
    struct ShuntData
    {
      using RealT = real_type;
      using IdxT  = index_type;

      std::string id;
      IdxT        bus{};
      RealT       G{};
      RealT       B{};
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

  } // namespace OPF
} // namespace GridKit
