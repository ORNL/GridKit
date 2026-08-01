#pragma once

#include <optional>
#include <string>

namespace GridKit
{
  namespace OPF
  {
    /// Immutable bus topology and parameter data.
    template <typename real_type, typename index_type>
    struct BusData
    {
      using RealT = real_type;
      using IdxT  = index_type;

      std::string          id;
      IdxT                 number{};
      RealT                kv{};
      std::optional<RealT> vmin;
      std::optional<RealT> vmax;
    };

  } // namespace OPF
} // namespace GridKit
