#pragma once

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/BusToSignalAdapter/BusToSignalAdapter.hpp>
#include <GridKit/Model/PhasorDynamics/BusToSignalAdapter/BusToSignalAdapterData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeData.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class BusToSignalAdapterTests
    {
    public:
      using AdapterT   = PhasorDynamics::BusToSignalAdapter<ScalarT, IdxT>;
      using RealT      = typename AdapterT::RealT;
      using BusT       = PhasorDynamics::Bus<ScalarT, IdxT>;
      using ComponentT = PhasorDynamics::Component<ScalarT, IdxT>;
      using SignalT    = PhasorDynamics::SignalNode<ScalarT, IdxT>;

      BusToSignalAdapterTests()  = default;
      ~BusToSignalAdapterTests() = default;

      TestOutcome constructor()
      {
        TestStatus success = true;

        auto* bus = new BusT(1.0, 0.0);

        ComponentT* adapter = new AdapterT(bus);

        success *= (adapter != nullptr);

        if (adapter)
        {
          delete adapter;
        }
        delete bus;

        return success.report(__func__);
      }

      TestOutcome residual()
      {
        TestStatus success = true;

        ScalarT Vr{10.0}; ///< Bus real voltage
        ScalarT Vi{20.0}; ///< Bus imaginary voltage

        const ScalarT Ir_ref{-5.0}; ///< Solution real current
        const ScalarT Ii_ref{0.0};  ///< Solution imaginary current
        ScalarT       Ir{Ir_ref};
        ScalarT       Ii{Ii_ref};

        auto bus = BusT(Vr, Vi);
        bus.allocate();
        bus.initialize();

        auto vr_sig = SignalT({.name = "vr", .signal_id = 0});
        auto vi_sig = SignalT({.name = "vi", .signal_id = 1});
        auto ir_sig = SignalT({.name = "ir", .signal_id = 2});
        auto ii_sig = SignalT({.name = "ii", .signal_id = 3});

        IdxT ir_index{INVALID_INDEX<IdxT>};
        IdxT ii_index{INVALID_INDEX<IdxT>};
        ir_sig.link(&Ir, &ir_index);
        ii_sig.link(&Ii, &ii_index);

        using namespace GridKit::PhasorDynamics;
        using SignalIn  = BusToSignalAdapterSignalInputs;
        using SignalOut = BusToSignalAdapterSignalOutputs;

        auto adapter = AdapterT(&bus);
        adapter.getPorts().out.template port<SignalOut::vr>().connect(&vr_sig);
        adapter.getPorts().out.template port<SignalOut::vi>().connect(&vi_sig);
        adapter.getPorts().in.template port<SignalIn::ir>().connect(&ir_sig);
        adapter.getPorts().in.template port<SignalIn::ii>().connect(&ii_sig);
        adapter.allocate();
        success *= (adapter.verify() == 0);
        success *= (vr_sig.read() == Vr);
        success *= (vi_sig.read() == Vi);

        adapter.evaluateResidual();
        success *= (bus.Ir() == Ir_ref);
        success *= (bus.Ii() == Ii_ref);

        return success.report(__func__);
      }
    };

  } // namespace Testing
} // namespace GridKit
