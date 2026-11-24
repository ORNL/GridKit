#pragma once
#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/ESST4B/Esst4bData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalPorts.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit::PhasorDynamics::Exciter
{
  template <class ScalarT, class IdxT>
  class Esst4b final : public Component<ScalarT, IdxT>
  {
    using RealT    = typename Component<ScalarT, IdxT>::RealT;
    using DataT    = Esst4bData<RealT, IdxT>;
    using Input    = typename DataT::SignalInputs;
    using Output   = typename DataT::SignalOutputs;
    using MonitorT = GridKit::Model::VariableMonitor<Esst4b, Esst4bData>;
    using BusT     = BusBase<ScalarT, IdxT>;

    enum : size_t
    {
      VC,
      XVR,
      VA,
      XVM,
      EFD,
      COUNT
    };

  public:
    explicit Esst4b(BusT* bus, const DataT& data);
    ~Esst4b() override;
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
    BusT*                       bus_;
    RealT                       vref_set_{};
    RealT                       voel_set_{};
    RealT                       kp_real_{}, kp_imag_{};
    RealT                       compensated_{};

    template <class T>
    [[gnu::always_inline]] T bridgeVoltage(T vr, T vi, T ir, T ii, T ifd) const;
    RealT                    mva_{};
    RealT                    Tr_{};
    RealT                    Kpr_{};
    RealT                    Kir_{};
    RealT                    Vrmax_{};
    RealT                    Vrmin_{};
    RealT                    Ta_{};
    RealT                    Kpm_{};
    RealT                    Kim_{};
    RealT                    Vmmax_{};
    RealT                    Vmmin_{};
    RealT                    Kg_{};
    RealT                    Kp_{};
    RealT                    Ki_{};
    RealT                    Vbmax_{};
    RealT                    Kc_{};
    RealT                    Xl_{};
    RealT                    ThetaP_{};
    RealT                    Vgmax_{};
  };
} // namespace GridKit::PhasorDynamics::Exciter
