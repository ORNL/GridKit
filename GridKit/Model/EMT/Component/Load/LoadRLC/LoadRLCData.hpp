/**
 * @file LoadRLCData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT parallel RLC loads.
 */

#pragma once

#include <cstddef>
#include <optional>
#include <tuple>
#include <vector>

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Terminals of a `LoadRLC`, each of size N
    enum class LoadRLCTerminals : size_t
    {
      bus, ///< Bus of the load
    };

    /// Parameters of a `LoadRLC`
    enum class LoadRLCParameters : size_t
    {
      N, ///< \f$N\f$ Number of phases
      R, ///< \f$R_n\f$ Resistance by phase [ohm]
      L, ///< \f$L_n\f$ Inductance by phase [H]
      C, ///< \f$C_n\f$ Capacitance by phase [F]
    };

    /// Data of a `LoadRLC`
    template <typename RealT, typename IdxT>
    struct LoadRLCData
    {
      using Parameters  = LoadRLCParameters;
      using PhaseValues = std::vector<std::optional<RealT>>; ///< No value for an absent branch

      IdxT                       N{0};
      std::optional<PhaseValues> R;
      std::optional<PhaseValues> L;
      std::optional<PhaseValues> C;

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::N, &LoadRLCData::N},
                                             Parameter{Parameters::R, &LoadRLCData::R},
                                             Parameter{Parameters::L, &LoadRLCData::L},
                                             Parameter{Parameters::C, &LoadRLCData::C}};
    };
  } // namespace EMT
} // namespace GridKit
