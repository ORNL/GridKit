/**
 * @file SwitchData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT switches.
 */

#pragma once

#include <cstddef>
#include <tuple>

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Ports of a `Switch`, each of size N
    enum class SwitchPorts : size_t
    {
      bus1, ///< Terminal 1
      bus2, ///< Terminal 2
    };

    /// Parameters of a `Switch`
    enum class SwitchParameters : size_t
    {
      N,    ///< \f$N\f$ Number of phases
      open, ///< \f$\mathrm{open}\f$ Initial command, `true` open and `false` closed
    };

    /// Data of a `Switch`
    template <typename RealT, typename IdxT>
    struct SwitchData
    {
      using Parameters = SwitchParameters;

      IdxT N{0};
      bool open{false};

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::N, &SwitchData::N}, Parameter{Parameters::open, &SwitchData::open}};
    };
  } // namespace EMT
} // namespace GridKit
