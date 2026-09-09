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
      I,          ///< Nominal phase RMS current for absolute tolerances [A]
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
      v1a, ///< Terminal 1 phase-a voltage
      v1b, ///< Terminal 1 phase-b voltage
      v1c, ///< Terminal 1 phase-c voltage
      v2a, ///< Terminal 2 phase-a voltage
      v2b, ///< Terminal 2 phase-b voltage
      v2c, ///< Terminal 2 phase-c voltage
      SIZE
    };

    /// Outputs supported by a lumped line
    enum class LineLumpedOutputs : size_t
    {
      i12a,
      i12b,
      i12c,
      i21a,
      i21b,
      i21c,
      SIZE
    };

    /// Variables able to be monitored for a lumped line
    enum class LineLumpedMonitorableVariables
    {
      i12a,
      i12b,
      i12c
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
