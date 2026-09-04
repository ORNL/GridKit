/**
 * @file VoltageSourceData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT voltage sources
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
    /// Initial parameters for a voltage source
    enum class VoltageSourceParameters
    {
      N,     ///< Number of phases
      E,     ///< Source voltage magnitudes, RMS
      phi,   ///< Source phase offsets
      omega, ///< Source angular frequency
      Rs,    ///< Series resistance matrix
      Ls,    ///< Series inductance matrix
    };

    /// Inputs supported by a voltage source
    enum class VoltageSourceInputs : size_t
    {
      bus, ///< Component ID of the connected bus
      SIZE
    };

    /// Outputs supported by a voltage source
    enum class VoltageSourceOutputs : size_t
    {
      SIZE
    };

    /// Variables able to be monitored for a voltage source
    enum class VoltageSourceMonitorableVariables
    {
      ea,
      eb,
      ec,
      ia,
      ib,
      ic
    };

    /**
     * @brief Contains modeling data for a voltage source
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     */
    template <typename real_type, typename index_type>
    struct VoltageSourceData : public ComponentData<real_type,
                                                    index_type,
                                                    VoltageSourceParameters,
                                                    VoltageSourceInputs,
                                                    VoltageSourceOutputs,
                                                    VoltageSourceMonitorableVariables>
    {
      VoltageSourceData() = default;

      using Parameters           = VoltageSourceParameters;
      using Inputs               = VoltageSourceInputs;
      using Outputs              = VoltageSourceOutputs;
      using MonitorableVariables = VoltageSourceMonitorableVariables;

      /// Rational source admittance submodel, replacing the series matrices
      /// when present
      std::optional<VectorFitData<real_type, index_type>> Y;
    };
  } // namespace EMT
} // namespace GridKit
