#pragma once

#include <GridKit/Model/PhasorDynamics/ModelData.hpp>
#include <GridKit/Model/PhasorDynamics/PortGroup.hpp>
#include <GridKit/Model/PhasorDynamics/SignalIn.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeSet.hpp>
#include <GridKit/Model/PhasorDynamics/SignalOut.hpp>

namespace GridKit::PhasorDynamics
{
  /// Set of input and output signal ports for a component.
  template <typename scalar_type, ModelData model_data_type>
  struct SignalPorts
  {
    using ScalarT        = scalar_type;
    using ModelDataT     = model_data_type;
    using IdxT           = ModelDataT::IdxT;
    using SignalInEnumT  = ModelDataT::SignalInputs;
    using SignalOutEnumT = ModelDataT::SignalOutputs;
    using SignalInT      = SignalIn<ScalarT, IdxT>;
    using SignalOutT     = SignalOut<ScalarT, IdxT>;
    using SignalNodeSetT = SignalNodeSet<ScalarT, IdxT>;

    PortGroup<SignalInT, SignalInEnumT>   in;
    PortGroup<SignalOutT, SignalOutEnumT> out;

    SignalPorts() = default;

    /// Connect this component's ports using the signal IDs in its model data.
    void connect(const ModelDataT& data, SignalNodeSetT& signal_nodes)
    {
      for (auto& [variable, id] : data.signal_inputs)
      {
        in[variable].connect(signal_nodes[id]);
      }
      for (auto& [variable, id] : data.signal_outputs)
      {
        out[variable].connect(signal_nodes[id]);
      }
    }
  };
} // namespace GridKit::PhasorDynamics
