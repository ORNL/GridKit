#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/SignalDelay/SignalDelayData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalHistory.hpp>
#include <GridKit/Model/PhasorDynamics/SignalPorts.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit::PhasorDynamics
{
  /// Algebraic delayed signal, created by the system composer.
  template <class ScalarT, class IdxT>
  class SignalDelay : public Component<ScalarT, IdxT>
  {
  public:
    using RealT        = typename Component<ScalarT, IdxT>::RealT;
    using MonitorT     = Model::VariableMonitor<SignalDelay, SignalDelayData>;
    using ModelDataT   = SignalDelayData<RealT, IdxT>;
    using SignalPortsT = SignalPorts<ScalarT, ModelDataT>;

    SignalDelay(const ModelDataT& data, SignalHistory<RealT>& history);
    ~SignalDelay() override;
    int  setGridKitComponentID(IdxT id) override;
    int  allocate() override;
    int  verify() const override;
    int  initialize() override;
    int  tagDifferentiable() override;
    int  setAbsoluteTolerance(RealT tolerance) override;
    int  evaluateResidual() override;
    int  evaluateJacobian() override;
    void linkOutput();

    SignalPortsT& getPorts()
    {
      return ports_;
    }

    const Model::VariableMonitorBase* getMonitor() const override;

  private:
    RealT                     delay_;
    SignalHistory<RealT>&     history_;
    SignalPortsT              ports_;
    std::unique_ptr<MonitorT> monitor_;
  };
} // namespace GridKit::PhasorDynamics
