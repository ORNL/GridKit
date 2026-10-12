#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit::PhasorDynamics
{
  enum class SignalDelayParameters : size_t
  {
    tau
  };
  enum class SignalDelayBuses : size_t
  {
  };
  enum class SignalDelayInputs : size_t
  {
    u
  };
  enum class SignalDelayOutputs : size_t
  {
    z
  };
  enum class SignalDelayMonitors : size_t
  {
    z
  };

  template <class RealT, class IdxT>
  using SignalDelayData = ComponentData<RealT, IdxT, SignalDelayParameters, SignalDelayBuses, SignalDelayInputs, SignalDelayOutputs, SignalDelayMonitors>;
} // namespace GridKit::PhasorDynamics
