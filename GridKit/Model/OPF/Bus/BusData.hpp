#pragma once

#include <optional>
#include <string>

namespace GridKit
{
  namespace OPF
  {
    /// Supported bus classes in the OPF input format.
    enum class BusClass
    {
      BUS,
      SLACK
    };

    /// Immutable bus topology and parameter data.
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

  } // namespace OPF
} // namespace GridKit
