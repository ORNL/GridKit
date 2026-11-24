#pragma once
#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit::PhasorDynamics::Governor
{
  enum class Ieeeg1Parameters : size_t
  {
    mva,
    Trate,
    K,
    T1,
    T2,
    T3,
    Uo,
    Uc,
    Pmax,
    Pmin,
    T4,
    K1,
    K2,
    T5,
    K3,
    K4,
    T6,
    K5,
    K6,
    T7,
    K7,
    K8,
    db1,
    db2,
    Eps,
    Gv1,
    Gv2,
    Gv3,
    Gv4,
    Gv5,
    Gv6,
    Pgv1,
    Pgv2,
    Pgv3,
    Pgv4,
    Pgv5,
    Pgv6
  };

  enum class Ieeeg1SignalInputs : size_t
  {
    speed,
    pref,
    paux
  };

  enum class Ieeeg1SignalOutputs : size_t
  {
    pmech_hp,
    pmech_lp
  };

  enum class Ieeeg1MonitorableVariables : size_t
  {
    pmech_hp,
    pmech_lp,
    pgv
  };

  template <class RealT, class IdxT>
  using Ieeeg1Data = ComponentData<RealT, IdxT, Ieeeg1Parameters, Utilities::EmptyEnum, Ieeeg1SignalInputs, Ieeeg1SignalOutputs, Ieeeg1MonitorableVariables, Utilities::EmptyEnum>;
} // namespace GridKit::PhasorDynamics::Governor
