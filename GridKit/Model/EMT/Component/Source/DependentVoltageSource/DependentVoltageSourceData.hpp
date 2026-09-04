/**
 * @file DependentVoltageSourceData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT dependent voltage sources
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
    /// Initial parameters for a dependent voltage source
    enum class DependentVoltageSourceParameters
    {
      N,  ///< Number of phases
      Rs, ///< Series resistance matrix
      Ls, ///< Series inductance matrix
    };

    /// Inputs supported by a dependent voltage source
    enum class DependentVoltageSourceInputs : size_t
    {
      bus, ///< Component ID of the connected bus
      ea,  ///< Phase-a source-voltage signal ID
      eb,  ///< Phase-b source-voltage signal ID
      ec,  ///< Phase-c source-voltage signal ID
      SIZE
    };

    /// Outputs supported by a dependent voltage source
    enum class DependentVoltageSourceOutputs : size_t
    {
      SIZE
    };

    /// Variables able to be monitored for a dependent voltage source
    enum class DependentVoltageSourceMonitorableVariables
    {
      ea,
      eb,
      ec,
      ia,
      ib,
      ic
    };

    /**
     * @brief Contains modeling data for a dependent voltage source
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     */
    template <typename real_type, typename index_type>
    struct DependentVoltageSourceData : public ComponentData<real_type,
                                                             index_type,
                                                             DependentVoltageSourceParameters,
                                                             DependentVoltageSourceInputs,
                                                             DependentVoltageSourceOutputs,
                                                             DependentVoltageSourceMonitorableVariables>
    {
      DependentVoltageSourceData() = default;

      using Parameters           = DependentVoltageSourceParameters;
      using Inputs               = DependentVoltageSourceInputs;
      using Outputs              = DependentVoltageSourceOutputs;
      using MonitorableVariables = DependentVoltageSourceMonitorableVariables;

      /// Rational source admittance submodel, replacing the series matrices
      /// when present
      std::optional<VectorFitData<real_type, index_type>> Y;
    };
  } // namespace EMT
} // namespace GridKit
