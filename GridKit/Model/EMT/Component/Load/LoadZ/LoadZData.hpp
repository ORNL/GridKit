/**
 * @file LoadZData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT impedance loads
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
    /// Initial parameters for an impedance load
    enum class LoadZParameters
    {
      N, ///< Number of phases
      R, ///< Load resistance matrix
      L, ///< Load inductance matrix
    };

    /// Inputs supported by an impedance load
    enum class LoadZInputs : size_t
    {
      bus, ///< Component ID of the connected bus
      SIZE
    };

    /// Outputs supported by an impedance load
    enum class LoadZOutputs : size_t
    {
      SIZE
    };

    /// Variables able to be monitored for an impedance load
    enum class LoadZMonitorableVariables
    {
      ia,
      ib,
      ic
    };

    /**
     * @brief Contains modeling data for an impedance load
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     */
    template <typename real_type, typename index_type>
    struct LoadZData : public ComponentData<real_type,
                                            index_type,
                                            LoadZParameters,
                                            LoadZInputs,
                                            LoadZOutputs,
                                            LoadZMonitorableVariables>
    {
      LoadZData() = default;

      using Parameters           = LoadZParameters;
      using Inputs               = LoadZInputs;
      using Outputs              = LoadZOutputs;
      using MonitorableVariables = LoadZMonitorableVariables;

      /// Rational impedance submodel, replacing the resistance and
      /// inductance matrices when present
      std::optional<VectorFitData<real_type, index_type>> Z;
    };
  } // namespace EMT
} // namespace GridKit
