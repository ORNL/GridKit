#pragma once

#include <optional>
#include <string>

namespace GridKit
{
  namespace OPF
  {
    /// Immutable branch topology and parameter data.
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

  } // namespace OPF
} // namespace GridKit
