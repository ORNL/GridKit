#pragma once

#include <string>

namespace GridKit
{
  namespace OPF
  {
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

  } // namespace OPF
} // namespace GridKit
