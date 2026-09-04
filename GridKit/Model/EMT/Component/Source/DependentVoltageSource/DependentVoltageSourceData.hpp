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

    /// Buses for a dependent voltage source
    enum class DependentVoltageSourceBuses : size_t
    {
      bus, ///< Unique ID of the bus to which the source is connected
      SIZE
    };

    /// Signal inputs supported for a dependent voltage source
    enum class DependentVoltageSourceSignalInputs : size_t
    {
      ea, ///< Phase-a source voltage
      eb, ///< Phase-b source voltage
      ec, ///< Phase-c source voltage
      SIZE
    };

    /// Signal outputs supported for a dependent voltage source
    enum class DependentVoltageSourceSignalOutputs : size_t
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
                                                             DependentVoltageSourceBuses,
                                                             DependentVoltageSourceSignalInputs,
                                                             DependentVoltageSourceSignalOutputs,
                                                             DependentVoltageSourceMonitorableVariables>
    {
      DependentVoltageSourceData() = default;

      using Parameters           = DependentVoltageSourceParameters;
      using Buses                = DependentVoltageSourceBuses;
      using SignalInputs         = DependentVoltageSourceSignalInputs;
      using SignalOutputs        = DependentVoltageSourceSignalOutputs;
      using MonitorableVariables = DependentVoltageSourceMonitorableVariables;

      /// Rational source admittance submodel, replacing the series matrices
      /// when present
      std::optional<VectorFitData<real_type, index_type>> Y;
    };
  } // namespace EMT
} // namespace GridKit
