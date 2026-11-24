#pragma once
#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/ESST1A/Esst1aData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalPorts.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit::PhasorDynamics::Exciter
{
  template <class ScalarT, class IdxT>
  class Esst1a final : public Component<ScalarT, IdxT>
  {
    using RealT    = typename Component<ScalarT, IdxT>::RealT;
    using DataT    = Esst1aData<RealT, IdxT>;
    using Input    = typename DataT::SignalInputs;
    using Output   = typename DataT::SignalOutputs;
    using MonitorT = GridKit::Model::VariableMonitor<Esst1a, Esst1aData>;
    using BusT     = BusBase<ScalarT, IdxT>;

    enum : size_t
    {
      VC,
      LL,
      LL1,
      VA,
      WASH,
      EFD,
      REG,
      DRIVE,
      STAGE,
      DRIVER,
      COUNT
    };

  public:
    explicit Esst1a(BusT* bus, const DataT& data);
    ~Esst1a() override;
    int allocate() override;
    int initialize() override;
    int evaluateResidual() override;
    int setSignalIndices() override;

    int setGridKitComponentID(IdxT id) override
    {
      this->gridkit_component_id_ = id;
      return 0;
    }

    int setAbsoluteTolerance(RealT tol) override
    {
      this->abs_tol_.setToConst(ScalarT{tol});
      return 0;
    }

    SignalPorts<ScalarT, DataT>& getPorts()
    {
      return ports_;
    }

    const GridKit::Model::VariableMonitorBase* getMonitor() const override
    {
      return monitor_.get();
    }

    GridKit::Model::ConfigurationChecks verify() const override;
    int                                 tagDifferentiable() override;
    int                                 evaluateJacobian() override;
    [[gnu::always_inline]] int          evaluateInternalResidual(const ScalarT* y, const ScalarT* yp, const ScalarT* wb, const ScalarT* w, ScalarT* f);

  private:
    SignalPorts<ScalarT, DataT> ports_;
    std::unique_ptr<MonitorT>   monitor_;
    RealT                       voel_set_{};
    RealT                       gate_input_{}, gate_output_{};
    RealT                       vref_set_{};
    BusT*                       bus_;

    RealT Tr_{};
    RealT Vimax_{};
    RealT Vimin_{};
    RealT Tc_{};
    RealT Tb_{};
    RealT Ka_{};
    RealT Ta_{};
    RealT Vrmax_{};
    RealT Vrmin_{};
    RealT Kc_{};
    RealT Kf_{};
    RealT Tf_{};
    IdxT  UEL_{1};
    IdxT  VOS_{1};
    RealT Tc1_{};
    RealT Tb1_{};
    RealT Vamax_{};
    RealT Vamin_{};
    RealT Klr_{};
    RealT Ilr_{};

    RealT washout_gain_, lead_gain_;
    RealT second_lead_gain_;
  };
} // namespace GridKit::PhasorDynamics::Exciter
