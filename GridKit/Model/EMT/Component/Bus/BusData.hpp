/**
 * @file BusData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT buses
 *
 */
#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>
#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFitData.hpp>

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
      ia,
      ib,
      ic,
      SIZE
    };

    /// Outputs supported by a bus
    enum class BusOutputs : size_t
    {
      va,
      vb,
      vc,
      SIZE
    };

    /// Variables able to be monitored for a bus
    enum class BusMonitorableVariables
    {
      va,
      vb,
      vc,
      i_sha,
      i_shb,
      i_shc
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

      /// Named bus shunts, each with an independent realization.
      std::map<std::string, VectorFitData<RealT, IdxT>> shunts;
    };
  } // namespace EMT
} // namespace GridKit
