/**
 * @file TransformerData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT transformers
 *
 */
#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Initial parameters for a transformer
    enum class TransformerParameters
    {
      N,     ///< Number of phases
      S,     ///< Rated three-phase apparent power
      V1,    ///< Rated line-to-line RMS voltage of winding 1
      V2,    ///< Rated line-to-line RMS voltage of winding 2
      f,     ///< Rated frequency
      P1,    ///< Terminal 1 connection map
      P2,    ///< Terminal 2 connection map
      tap,   ///< Off-nominal ratio on winding 1
      R,     ///< Short-circuit resistance
      X,     ///< Short-circuit reactance
      I0,    ///< No-load current at rated voltage
      P0,    ///< No-load loss at rated voltage
      knee,  ///< Knee flux linkage
      Lsat,  ///< Terminal saturation inductance
      split, ///< Winding 1 share of the magnetizing branch
    };

    /// Inputs supported by a transformer
    enum class TransformerInputs : size_t
    {
      v1a, ///< Terminal 1 phase-a voltage
      v1b, ///< Terminal 1 phase-b voltage
      v1c, ///< Terminal 1 phase-c voltage
      v2a, ///< Terminal 2 phase-a voltage
      v2b, ///< Terminal 2 phase-b voltage
      v2c, ///< Terminal 2 phase-c voltage
      SIZE
    };

    /// Outputs supported by a transformer
    enum class TransformerOutputs : size_t
    {
      i1a, ///< Terminal 1 phase-a current injection [A]
      i1b, ///< Terminal 1 phase-b current injection [A]
      i1c, ///< Terminal 1 phase-c current injection [A]
      i2a, ///< Terminal 2 phase-a current injection [A]
      i2b, ///< Terminal 2 phase-b current injection [A]
      i2c, ///< Terminal 2 phase-c current injection [A]
      SIZE
    };

    /// Variables able to be monitored for a transformer
    enum class TransformerMonitorableVariables
    {
      i12a,
      i12b,
      i12c,
      psi1a,
      psi1b,
      psi1c,
      psi2a,
      psi2b,
      psi2c,
      i1a,
      i1b,
      i1c,
      i2a,
      i2b,
      i2c
    };

    /**
     * @brief Contains modeling data for a transformer
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     */
    template <typename real_type, typename index_type>
    struct TransformerData : public ComponentData<real_type,
                                                  index_type,
                                                  TransformerParameters,
                                                  TransformerInputs,
                                                  TransformerOutputs,
                                                  TransformerMonitorableVariables>
    {
      TransformerData() = default;

      using Parameters           = TransformerParameters;
      using Inputs               = TransformerInputs;
      using Outputs              = TransformerOutputs;
      using MonitorableVariables = TransformerMonitorableVariables;
    };
  } // namespace EMT
} // namespace GridKit
