#pragma once

#include <optional>
#include <string>

namespace GridKit
{
  namespace OPF
  {
    /// Immutable generator topology, limits, and cost data.
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

  } // namespace OPF
} // namespace GridKit
