#pragma once
#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/Controller/REECA/ReecaData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalPorts.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit::PhasorDynamics::Controller
{
  template <class ScalarT, class IdxT>
  class Reeca final : public Component<ScalarT, IdxT>
  {
    using RealT    = typename Component<ScalarT, IdxT>::RealT;
    using DataT    = ReecaData<RealT, IdxT>;
    using Input    = typename DataT::SignalInputs;
    using Output   = typename DataT::SignalOutputs;
    using MonitorT = GridKit::Model::VariableMonitor<Reeca, ReecaData>;
    using BusT     = BusBase<ScalarT, IdxT>;

    enum : size_t
    {
      VTERM,
      VMEAS,
      PMEAS,
      XPIQ,
      XPIV,
      QV,
      PREF_RATE,
      PORD,
      IQV,
      IPMAX,
      IQCMD,
      IPCMD,
      IQMAX,
      ILCAP,
      VSAFE,
      EQ,
      EPIV,
      COUNT
    };

  public:
    explicit Reeca(BusT* bus, const DataT& data);
    ~Reeca() override;
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

    IdxT sizeRoots() override
    {
      return 4;
    }

    int  evaluateRootsAt(RealT t, const ScalarT* y, RealT* roots) override;
    bool applyRootEvents(RealT t, std::span<const int> crossed) override;
    bool synchronizeRootEvents(RealT t, const ScalarT* previous) override;

  private:
    SignalPorts<ScalarT, DataT> ports_;
    std::unique_ptr<MonitorT>   monitor_;
    RealT                       pe_set_{};
    RealT                       pfaref_set_{};
    RealT                       pref_set_{};
    RealT                       qext_set_{};
    RealT                       qgen_set_{};
    BusT*                       bus_;

    template <class T>
    [[gnu::always_inline]] T currentLimit(T voltage, const std::array<RealT, 4>& v, const std::array<RealT, 4>& i, const std::array<RealT, 3>& slope) const;
    RealT                    mva_{ZERO<RealT>};
    bool                     PfFlag_{};
    bool                     VFlag_{};
    bool                     QFlag_{};
    bool                     Pqflag_{};
    bool                     Pflag_{};
    RealT                    Trv_{0.02};
    RealT                    Tp_{ZERO<RealT>};
    RealT                    Vref0_{ZERO<RealT>};
    RealT                    Vref1_{};
    RealT                    Vdip_{0.85};
    RealT                    Vup_{1.15};
    RealT                    dbd1_{ZERO<RealT>};
    RealT                    dbd2_{ZERO<RealT>};
    RealT                    kqv_{5.0};
    RealT                    Iql1_{-1.1};
    RealT                    Iqh1_{1.1};
    RealT                    Iqfrz_{};
    RealT                    Thld_{};
    RealT                    Thld2_{};
    RealT                    Qmax_{0.436};
    RealT                    Qmin_{-0.436};
    RealT                    Kqp_{ZERO<RealT>};
    RealT                    Kqi_{0.1};
    RealT                    Vmax_{1.1};
    RealT                    Vmin_{0.9};
    RealT                    Kvp_{18.0};
    RealT                    Kvi_{5.0};
    RealT                    Tiq_{0.02};
    RealT                    Tpord_{0.02};
    RealT                    dPmax_{99.0};
    RealT                    dPmin_{-99.0};
    RealT                    Pmax_{ONE<RealT>};
    RealT                    Pmin_{ZERO<RealT>};
    RealT                    Imax_{1.3};
    std::array<RealT, 4>     Vq_{};
    std::array<RealT, 4>     Iq_{};
    std::array<RealT, 4>     Vp_{};
    std::array<RealT, 4>     Ip_{};

    // Slopes of the three voltage-dependent current limit segments.
    std::array<RealT, 3> slope_q_{}, slope_p_{};
    void                 voltageDip(bool dip, RealT t, RealT ipmax);
    RealT                reference_voltage_{0};
    RealT                dip_{0};
    RealT                hold_p_{0}, hold_q_{0};
    RealT                release_p_{0}, release_q_{0};
    RealT                held_ipmax_{0};
  };
} // namespace GridKit::PhasorDynamics::Controller
