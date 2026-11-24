#pragma once
#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/GGOV1/Ggov1Data.hpp>
#include <GridKit/Model/PhasorDynamics/SignalPorts.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit::PhasorDynamics::Governor
{
  template <class ScalarT, class IdxT>
  class Ggov1 final : public Component<ScalarT, IdxT>
  {
    using RealT    = typename Component<ScalarT, IdxT>::RealT;
    using DataT    = Ggov1Data<RealT, IdxT>;
    using Input    = typename DataT::SignalInputs;
    using Output   = typename DataT::SignalOutputs;
    using MonitorT = GridKit::Model::VariableMonitor<Ggov1, Ggov1Data>;

    enum : size_t
    {
      PE,
      XD,
      XI,
      VALVE,
      TURB,
      TEMP,
      LOAD,
      XL,
      XMW,
      SPEED,
      DEMAND,
      FSR,
      PM,
      ERROR,
      MEMORY, // Fuel-command memory deviation from initial_fsr_.
      FSRN,
      FSRT,
      COUNT
    };

  public:
    explicit Ggov1(const DataT& data);
    ~Ggov1() override;
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
    RealT                       pmwset_set_{};
    RealT                       pref_set_{};
    RealT                       initial_fsr_{};

    RealT                    mva_{};
    RealT                    Trate_{};
    std::make_signed_t<IdxT> Rselect_{1};
    bool                     Flag_{};
    RealT                    R_{};
    RealT                    Tpelec_{};
    RealT                    Maxerr_{};
    RealT                    Minerr_{};
    RealT                    Kpgov_{};
    RealT                    Kigov_{};
    RealT                    Kdgov_{};
    RealT                    Tdgov_{};
    RealT                    Vmax_{};
    RealT                    Vmin_{};
    RealT                    Tact_{};
    RealT                    Kturb_{};
    RealT                    Wfnl_{};
    RealT                    Tb_{};
    RealT                    Tc_{};
    RealT                    Teng_{};
    RealT                    Tfload_{};
    RealT                    Kpload_{};
    RealT                    Kiload_{};
    RealT                    Ldref_{};
    RealT                    Dm_{};
    RealT                    Ropen_{};
    RealT                    Rclose_{};
    RealT                    Kimw_{};
    RealT                    Aset_{};
    RealT                    Ka_{};
    RealT                    Ta_{};
    RealT                    Db_{};
    RealT                    Tsa_{};
    RealT                    Tsb_{};
    RealT                    Tfsr_{1e-3};
    RealT                    Fsrtmax_{1};

    RealT governor_tracking_, load_tracking_, temperature_gain_, turbine_gain_, acceleration_gain_, derivative_gain_;
  };
} // namespace GridKit::PhasorDynamics::Governor
