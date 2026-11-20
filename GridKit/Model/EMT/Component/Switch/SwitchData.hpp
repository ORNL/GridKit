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
    /// Terminals of a `Switch`, each of size N
    enum class SwitchTerminals : size_t
    {
      bus1, ///< Terminal 1
      bus2, ///< Terminal 2
    };

    /// Inputs of a `Switch`
    enum class SwitchInputs : size_t
    {
      open, ///< \f$\mathrm{open}\f$ Command by phase, 1 open and 0 closed
    };

    /// Parameters of a `Switch`
    enum class SwitchParameters : size_t
    {
      N, ///< \f$N\f$ Number of phases
    };

    /// Data of a `Switch`
    template <typename RealT, typename IdxT>
    struct SwitchData
    {
      using Parameters = SwitchParameters;

      IdxT N{0};

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::N, &SwitchData::N}};
    };
  } // namespace EMT
} // namespace GridKit
