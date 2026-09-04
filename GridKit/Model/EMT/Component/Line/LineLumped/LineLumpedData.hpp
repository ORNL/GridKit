/**
 * @file LineLumpedData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT lumped lines
 *
 */
#pragma once

#include <optional>

#include <GridKit/Model/EMT/ComponentData.hpp>
#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFitData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Initial parameters for a lumped line
    enum class LineLumpedParameters
    {
      N,          ///< Number of phases
      K,          ///< Number of conductors
      conductors, ///< Conductor phase-index list
      dx,         ///< Line segment length
      Rp,         ///< Series resistance matrix per unit length
      Lp,         ///< Series inductance matrix per unit length
      Gp,         ///< Shunt conductance matrix per unit length
      Cp,         ///< Shunt capacitance matrix per unit length
    };

    /// Inputs supported by a lumped line
    enum class LineLumpedInputs : size_t
    {
      bus1, ///< Component ID of the terminal 1 bus
      bus2, ///< Component ID of the terminal 2 bus
      SIZE
    };

    /// Outputs supported by a lumped line
    enum class LineLumpedOutputs : size_t
    {
      SIZE
    };

    /// Variables able to be monitored for a lumped line
    enum class LineLumpedMonitorableVariables
    {
      i12a,
      i12b,
      i12c,
      i_sh1a,
      i_sh1b,
      i_sh1c,
      i_sh2a,
      i_sh2b,
      i_sh2c
    };

    /**
     * @brief Contains modeling data for a lumped line
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     */
    template <typename real_type, typename index_type>
    struct LineLumpedData : public ComponentData<real_type,
                                                 index_type,
                                                 LineLumpedParameters,
                                                 LineLumpedInputs,
                                                 LineLumpedOutputs,
                                                 LineLumpedMonitorableVariables>
    {
      LineLumpedData() = default;

      using Parameters           = LineLumpedParameters;
      using Inputs               = LineLumpedInputs;
      using Outputs              = LineLumpedOutputs;
      using MonitorableVariables = LineLumpedMonitorableVariables;

      /// Rational per-unit-length series impedance submodel, replacing the
      /// series matrices when present
      std::optional<VectorFitData<real_type, index_type>> Zp;

      /// Rational per-unit-length shunt admittance coefficient set,
      /// instantiated once per terminal, replacing the shunt matrices when
      /// present
      std::optional<VectorFitData<real_type, index_type>> Yp;
    };
  } // namespace EMT
} // namespace GridKit
