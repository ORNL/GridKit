/**
 * @file BusData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT buses
 *
 */
#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Initial parameters for a bus
    enum class BusParameters
    {
      N, ///< Number of phases
    };

    /// Inputs supported by a bus
    enum class BusInputs : size_t
    {
      SIZE
    };

    /// Outputs supported by a bus
    enum class BusOutputs : size_t
    {
      SIZE
    };

    /// Variables able to be monitored for a bus
    enum class BusMonitorableVariables
    {
      va,
      vb,
      vc
    };

    /**
     * @brief Contains modeling data for a bus
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     */
    template <typename real_type, typename index_type>
    struct BusData : public ComponentData<real_type,
                                          index_type,
                                          BusParameters,
                                          BusInputs,
                                          BusOutputs,
                                          BusMonitorableVariables>
    {
      BusData() = default;

      using Parameters           = BusParameters;
      using Inputs               = BusInputs;
      using Outputs              = BusOutputs;
      using MonitorableVariables = BusMonitorableVariables;

      using IdxT  = index_type;
      using RealT = real_type;

      RealT va0{0.0}; ///< Initial instantaneous phase a voltage
      RealT vb0{0.0}; ///< Initial instantaneous phase b voltage
      RealT vc0{0.0}; ///< Initial instantaneous phase c voltage
    };
  } // namespace EMT
} // namespace GridKit
