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

    /// Buses for an impedance load
    enum class LoadZBuses : size_t
    {
      bus, ///< Unique ID of the bus to which the load is connected
      SIZE
    };

    /// Signal inputs supported for an impedance load
    enum class LoadZSignalInputs : size_t
    {
      SIZE
    };

    /// Signal outputs supported for an impedance load
    enum class LoadZSignalOutputs : size_t
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
                                            LoadZBuses,
                                            LoadZSignalInputs,
                                            LoadZSignalOutputs,
                                            LoadZMonitorableVariables>
    {
      LoadZData() = default;

      using Parameters           = LoadZParameters;
      using Buses                = LoadZBuses;
      using SignalInputs         = LoadZSignalInputs;
      using SignalOutputs        = LoadZSignalOutputs;
      using MonitorableVariables = LoadZMonitorableVariables;

      /// Rational impedance submodel, replacing the resistance and
      /// inductance matrices when present
      std::optional<VectorFitData<real_type, index_type>> Z;
    };
  } // namespace EMT
} // namespace GridKit
