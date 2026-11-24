#pragma once
#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit::PhasorDynamics::Exciter
{
  enum class Esst4bParameters : size_t
  {
    mva,
    Tr,
    Kpr,
    Kir,
    Vrmax,
    Vrmin,
    Ta,
    Kpm,
    Kim,
    Vmmax,
    Vmmin,
    Kg,
    Kp,
    Ki,
    Vbmax,
    Kc,
    Xl,
    ThetaP,
    Vgmax
  };

  enum class Esst4bBuses : size_t
  {
    bus
  };

  enum class Esst4bSignalInputs : size_t
  {
    ifd,
    ir,
    ii,
    vref,
    vs,
    vuel,
    voel,
    vcomp
  };

  enum class Esst4bSignalOutputs : size_t
  {
    efd
  };

  enum class Esst4bMonitorableVariables : size_t
  {
    efd
  };

  template <class RealT, class IdxT>
  using Esst4bData = ComponentData<RealT, IdxT, Esst4bParameters, Esst4bBuses, Esst4bSignalInputs, Esst4bSignalOutputs, Esst4bMonitorableVariables, Utilities::EmptyEnum>;
} // namespace GridKit::PhasorDynamics::Exciter
