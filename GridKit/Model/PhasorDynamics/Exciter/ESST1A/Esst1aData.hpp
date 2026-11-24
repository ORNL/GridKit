#pragma once
#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit::PhasorDynamics::Exciter
{
  enum class Esst1aParameters : size_t
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
    Tf,
    UEL,
    VOS,
    Tc1,
    Tb1,
    Vamax,
    Vamin,
    Klr,
    Ilr
  };

  enum class Esst1aBuses : size_t
  {
    bus
  };

  enum class Esst1aSignalInputs : size_t
  {
    ifd,
    vref,
    vs,
    vuel,
    voel
  };

  enum class Esst1aSignalOutputs : size_t
  {
    efd
  };

  enum class Esst1aMonitorableVariables : size_t
  {
    efd
  };

  template <class RealT, class IdxT>
  using Esst1aData = ComponentData<RealT, IdxT, Esst1aParameters, Esst1aBuses, Esst1aSignalInputs, Esst1aSignalOutputs, Esst1aMonitorableVariables, Utilities::EmptyEnum>;
} // namespace GridKit::PhasorDynamics::Exciter
