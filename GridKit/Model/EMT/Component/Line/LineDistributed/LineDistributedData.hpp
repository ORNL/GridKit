#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>
#include <GridKit/Model/EMT/Operators/Shift/Propagation/PropagationData.hpp>

namespace GridKit
{
  namespace EMT
  {
    enum class LineDistributedParameters
    {
      N,          ///< Number of phases
      K,          ///< Number of conductors
      conductors, ///< Conductor phase-index list
    };
    enum class LineDistributedInputs : size_t
    {
      v1a, ///< Terminal 1 phase-a voltage
      v1b, ///< Terminal 1 phase-b voltage
      v1c, ///< Terminal 1 phase-c voltage
      v2a, ///< Terminal 2 phase-a voltage
      v2b, ///< Terminal 2 phase-b voltage
      v2c, ///< Terminal 2 phase-c voltage
      SIZE,
    };
    enum class LineDistributedOutputs : size_t
    {
      i_ref1a,
      i_ref1b,
      i_ref1c,
      i_ref2a,
      i_ref2b,
      i_ref2c,
      i_inc1a,
      i_inc1b,
      i_inc1c,
      i_inc2a,
      i_inc2b,
      i_inc2c,
      SIZE,
    };
    enum class LineDistributedMonitorableVariables
    {
      i_c1a,
      i_c1b,
      i_c1c,
      i_c2a,
      i_c2b,
      i_c2c,
      i_inc1a,
      i_inc1b,
      i_inc1c,
      i_inc2a,
      i_inc2b,
      i_inc2c,
      i_ref1a,
      i_ref1b,
      i_ref1c,
      i_ref2a,
      i_ref2b,
      i_ref2c,
    };

    template <typename real_type, typename index_type>
    struct LineDistributedData : ComponentData<real_type, index_type, LineDistributedParameters, LineDistributedInputs, LineDistributedOutputs, LineDistributedMonitorableVariables>
    {
      LineDistributedData() = default;

      using Parameters           = LineDistributedParameters;
      using Inputs               = LineDistributedInputs;
      using Outputs              = LineDistributedOutputs;
      using MonitorableVariables = LineDistributedMonitorableVariables;
      VectorFitData<real_type, index_type>   Yc;
      PropagationData<real_type, index_type> H;
    };
  } // namespace EMT
} // namespace GridKit
