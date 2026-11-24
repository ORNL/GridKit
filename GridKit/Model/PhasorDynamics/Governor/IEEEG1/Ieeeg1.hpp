#pragma once
#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/IEEEG1/Ieeeg1Data.hpp>
#include <GridKit/Model/PhasorDynamics/SignalPorts.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit::PhasorDynamics::Governor
{
  template <class ScalarT, class IdxT>
  class Ieeeg1 final : public Component<ScalarT, IdxT>
  {
    using RealT    = typename Component<ScalarT, IdxT>::RealT;
    using DataT    = Ieeeg1Data<RealT, IdxT>;
    using Input    = typename DataT::SignalInputs;
    using Output   = typename DataT::SignalOutputs;
    using MonitorT = GridKit::Model::VariableMonitor<Ieeeg1, Ieeeg1Data>;

    enum : size_t
    {
      SPEED,
      VALVE,
      BOWL,
      REHEAT,
      CROSSOVER,
      REHEAT2,
      PM,
      PM_LP,
      COUNT
    };

  public:
    explicit Ieeeg1(const DataT& data);
    ~Ieeeg1() override;
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
    RealT                       pref_set_{};

    template <class T>
    [[gnu::always_inline]] T valveFlow(T valve) const;
    RealT                    mva_{};
    RealT                    Trate_{};
    RealT                    K_{};
    RealT                    T1_{};
    RealT                    T2_{};
    RealT                    T3_{};
    RealT                    Uo_{};
    RealT                    Uc_{};
    RealT                    Pmax_{};
    RealT                    Pmin_{};
    RealT                    T4_{};
    RealT                    K1_{};
    RealT                    K2_{};
    RealT                    T5_{};
    RealT                    K3_{};
    RealT                    K4_{};
    RealT                    T6_{};
    RealT                    K5_{};
    RealT                    K6_{};
    RealT                    T7_{};
    RealT                    K7_{};
    RealT                    K8_{};
    RealT                    db1_{};
    RealT                    db2_{};
    RealT                    Eps_{};
    std::array<RealT, 6>     Gv_{};
    std::array<RealT, 6>     Pgv_{};

    RealT                lead_gain_{}, straight_valve_{};
    std::array<RealT, 5> slope_{};
  };
} // namespace GridKit::PhasorDynamics::Governor
