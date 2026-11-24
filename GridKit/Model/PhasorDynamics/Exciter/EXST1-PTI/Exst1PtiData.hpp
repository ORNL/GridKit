#pragma once
#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit::PhasorDynamics::Exciter
{
  enum class Exst1PtiParameters : size_t
  {
    Tr,
    Vimax,
    Vimin,
    Tc,
    Tb,
    Ka,
    Ta,
    Vrmax,
    Vrmin,
    Kc,
    Kf,
    Tf
  };

  enum class Exst1PtiBuses : size_t
  {
    bus
  };

  enum class Exst1PtiSignalInputs : size_t
  {
    ifd,
    vref,
    vs,
    vuel,
    voel
  };

  enum class Exst1PtiSignalOutputs : size_t
  {
    efd
  };

  enum class Exst1PtiMonitorableVariables : size_t
  {
    efd
  };

  template <class RealT, class IdxT>
  using Exst1PtiData = ComponentData<RealT, IdxT, Exst1PtiParameters, Exst1PtiBuses, Exst1PtiSignalInputs, Exst1PtiSignalOutputs, Exst1PtiMonitorableVariables, Utilities::EmptyEnum>;
} // namespace GridKit::PhasorDynamics::Exciter
