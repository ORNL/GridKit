#pragma once

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Connector/CoSim.hpp>
#include <GridKit/Model/PhasorDynamics/Connector/CoSimData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeData.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class CoSimTests
    {
    public:
      using CoSimT     = PhasorDynamics::Connector::CoSim<ScalarT, IdxT>;
      using RealT      = typename CoSimT::RealT;
      using BusT       = PhasorDynamics::Bus<ScalarT, IdxT>;
      using ComponentT = PhasorDynamics::Component<ScalarT, IdxT>;
      using SignalT    = PhasorDynamics::SignalNode<ScalarT, IdxT>;

      CoSimTests()  = default;
      ~CoSimTests() = default;

      TestOutcome constructor()
      {
        TestStatus success = true;

        auto* bus = new BusT(1.0, 0.0);

        ComponentT* cosim = new CoSimT(bus);

        success *= (cosim != nullptr);

        if (cosim)
        {
          delete cosim;
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

        IdxT ir_index;
        IdxT ii_index;
        ir_sig.set(&Ir, &ir_index);
        ii_sig.set(&Ii, &ii_index);

        using namespace GridKit::PhasorDynamics::Connector;
        constexpr auto VREAL = CoSimInternalVariables::VREAL;
        constexpr auto VIMAG = CoSimInternalVariables::VIMAG;
        constexpr auto IREAL = CoSimExternalVariables::IREAL;
        constexpr auto IIMAG = CoSimExternalVariables::IIMAG;

        auto cosim = CoSimT(&bus);
        cosim.getSignals().template assignSignalNode<VREAL>(&vr_sig);
        cosim.getSignals().template assignSignalNode<VIMAG>(&vi_sig);
        cosim.getSignals().template attachSignalNode<IREAL>(&ir_sig);
        cosim.getSignals().template attachSignalNode<IIMAG>(&ii_sig);
        cosim.allocate();
        success *= (cosim.verify() == 0);
        success *= (vr_sig.read() == Vr);
        success *= (vi_sig.read() == Vi);

        cosim.evaluateResidual();
        success *= (bus.Ir() == Ir_ref);
        success *= (bus.Ii() == Ii_ref);

        return success.report(__func__);
      }
    };

  } // namespace Testing
} // namespace GridKit
