#pragma once
#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit::PhasorDynamics::Controller
{
  enum class ReecaParameters : size_t
  {
    mva,
    PfFlag,
    VFlag,
    QFlag,
    Pqflag,
    Pflag,
    Trv,
    Tp,
    Vref0,
    Vref1,
    Vdip,
    Vup,
    dbd1,
    dbd2,
    kqv,
    Iql1,
    Iqh1,
    Iqfrz,
    Thld,
    Thld2,
    Qmax,
    Qmin,
    Kqp,
    Kqi,
    Vmax,
    Vmin,
    Kvp,
    Kvi,
    Tiq,
    Tpord,
    dPmax,
    dPmin,
    Pmax,
    Pmin,
    Imax,
    Vq1,
    Vq2,
    Vq3,
    Vq4,
    Iq1,
    Iq2,
    Iq3,
    Iq4,
    Vp1,
    Vp2,
    Vp3,
    Vp4,
    Ip1,
    Ip2,
    Ip3,
    Ip4
  };

  enum class ReecaBuses : size_t
  {
    bus
  };

  enum class ReecaSignalInputs : size_t
  {
    pe,
    qgen,
    qext,
    pfaref,
    pref,
    speed
  };

  enum class ReecaSignalOutputs : size_t
  {
    iqcmd,
    ipcmd
  };

  enum class ReecaMonitorableVariables : size_t
  {
    iqcmd,
    ipcmd,
    vt,
    ipmax,
    iqinj
  };

  template <class RealT, class IdxT>
  using ReecaData = ComponentData<RealT, IdxT, ReecaParameters, ReecaBuses, ReecaSignalInputs, ReecaSignalOutputs, ReecaMonitorableVariables, Utilities::EmptyEnum>;
} // namespace GridKit::PhasorDynamics::Controller
