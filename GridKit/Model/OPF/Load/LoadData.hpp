#pragma once

#include <optional>
#include <string>

namespace GridKit
{
  namespace OPF
  {
    /// Immutable load topology and limit data.
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

  } // namespace OPF
} // namespace GridKit
