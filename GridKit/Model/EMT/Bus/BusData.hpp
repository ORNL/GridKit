/**
 * @file BusData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT buses.
 */

#pragma once

#include <cstddef>
#include <tuple>

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Parameters of a `Bus`
    enum class BusParameters : size_t
    {
      N,  ///< \f$N\f$ Number of phases
      kv, ///< \f$V^\mathrm{base}\f$ Base voltage [kV]
    };

    /// Data of a `Bus`
    template <typename RealT, typename IdxT>
    struct BusData
    {
      using Parameters = BusParameters;

      IdxT  N{0};
      RealT kv{ZERO<RealT>};

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::N, &BusData::N}, Parameter{Parameters::kv, &BusData::kv}};
    };
  } // namespace EMT
} // namespace GridKit
