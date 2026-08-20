#pragma once

#include <concepts>

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit::PhasorDynamics
{
  /// Concept implemented by component-model data containers.
  template <typename T>
  concept ModelData = requires {
    typename T::RealT;
    typename T::IdxT;
    typename T::Parameters;
    typename T::Buses;
    typename T::SignalInputs;
    typename T::SignalOutputs;
    typename T::MonitorableVariables;
  } && std::derived_from<T, ComponentData<typename T::RealT, typename T::IdxT, typename T::Parameters, typename T::Buses, typename T::SignalInputs, typename T::SignalOutputs, typename T::MonitorableVariables>> && Utilities::SizedEnum<typename T::Parameters> && Utilities::SizedEnum<typename T::Buses> && Utilities::SizedEnum<typename T::SignalInputs> && Utilities::SizedEnum<typename T::SignalOutputs> && Utilities::SizedEnum<typename T::MonitorableVariables>;
} // namespace GridKit::PhasorDynamics
