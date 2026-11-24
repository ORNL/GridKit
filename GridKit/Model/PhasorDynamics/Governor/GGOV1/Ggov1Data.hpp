#pragma once
#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit::PhasorDynamics::Governor
{
  enum class Ggov1Parameters : size_t
  {
    mva,
    Trate,
    Rselect,
    Flag,
    R,
    Tpelec,
    Maxerr,
    Minerr,
    Kpgov,
    Kigov,
    Kdgov,
    Tdgov,
    Vmax,
    Vmin,
    Tact,
    Kturb,
    Wfnl,
    Tb,
    Tc,
    Teng,
    Tfload,
    Kpload,
    Kiload,
    Ldref,
    Dm,
    Ropen,
    Rclose,
    Kimw,
    Aset,
    Ka,
    Ta,
    Db,
    Tsa,
    Tsb,
    Tfsr,
    Fsrtmax
  };

  enum class Ggov1SignalInputs : size_t
  {
    speed,
    pelec,
    pref,
    paux,
    pmwset
  };

  enum class Ggov1SignalOutputs : size_t
  {
    pmech
  };

  enum class Ggov1MonitorableVariables : size_t
  {
    pmech,
    valve
  };

  template <class RealT, class IdxT>
  using Ggov1Data = ComponentData<RealT, IdxT, Ggov1Parameters, Utilities::EmptyEnum, Ggov1SignalInputs, Ggov1SignalOutputs, Ggov1MonitorableVariables, Utilities::EmptyEnum>;
} // namespace GridKit::PhasorDynamics::Governor
