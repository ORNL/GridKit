#pragma once

#include <iostream>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Constants.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusSignalVoltageIn/BusSignalVoltageIn.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusSignalVoltageIn/BusSignalVoltageInData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeData.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class BusSignalVoltageInTests
    {
    public:
      using BusT      = PhasorDynamics::BusSignalVoltageIn<ScalarT, IdxT>;
      using BusTypeT  = typename BusT::BusTypeT;
      using SignalT   = PhasorDynamics::SignalNode<ScalarT, IdxT>;
      using SignalIn  = PhasorDynamics::BusSignalVoltageInInputs;
      using SignalOut = PhasorDynamics::BusSignalVoltageInOutputs;

      BusSignalVoltageInTests()  = default;
      ~BusSignalVoltageInTests() = default;

      /// Constructor, allocation, and initialization checks
      TestOutcome constructor()
      {
        TestStatus success = true;

        const ScalarT Vr{0.93};
        const ScalarT Vi{-0.27};

        PhasorDynamics::BusBase<ScalarT, IdxT>* bus = nullptr;

        bus = new BusT();
        bus->allocate();
        bus->initialize();
        success *= isEqual(bus->Vr(), 0.0);
        success *= isEqual(bus->Vi(), 0.0);
        success *= (bus->size() == 0);
        success *= (bus->BusType() == BusTypeT::SIGNAL_VOLTAGE_IN);
        delete bus;

        bus      = new BusT(Vr, Vi);
        success *= isEqual(bus->Vr(), Vr);
        success *= isEqual(bus->Vi(), Vi);
        bus->allocate();
        bus->initialize();
        success *= isEqual(bus->Vr(), Vr);
        success *= isEqual(bus->Vi(), Vi);
        success *= isEqual(bus->Ir(), 0.0);
        success *= isEqual(bus->Ii(), 0.0);
        delete bus;

        bus = nullptr;

        return success.report(__func__);
      }

      /// Input ports set the bus voltage
      TestOutcome voltageInputs()
      {
        TestStatus success = true;

        const ScalarT Vr0{1.0};
        const ScalarT Vi0{0.0};
        ScalarT       Vr{0.93};  ///< Voltage on signal vr
        ScalarT       Vi{-0.27}; ///< Voltage on signal vi
        IdxT          vr_index{7};
        IdxT          vi_index{8};

        auto vr_node = SignalT({.name = "vr", .signal_id = 0});
        auto vi_node = SignalT({.name = "vi", .signal_id = 1});
        vr_node.link(&Vr, &vr_index);
        vi_node.link(&Vi, &vi_index);

        BusT bus(Vr0, Vi0);
        bus.getPorts().in.template port<SignalIn::vr>().connect(&vr_node);
        bus.getPorts().in.template port<SignalIn::vi>().connect(&vi_node);
        bus.allocate();
        bus.initialize();
        success *= (bus.verify() == 0);

        // Voltage is read straight from the signals, no evaluation needed
        success *= isEqual(bus.Vr(), Vr);
        success *= isEqual(bus.Vi(), Vi);
        success *= (&bus.Vr() == &Vr);
        success *= (&bus.Vi() == &Vi);

        // Voltage follows the signals
        Vr       = 1.17;
        Vi       = 0.41;
        success *= isEqual(bus.Vr(), 1.17);
        success *= isEqual(bus.Vi(), 0.41);
        bus.evaluateResidual();
        success *= isEqual(bus.Vr(), 1.17);
        success *= isEqual(bus.Vi(), 0.41);

        // Unconnected inputs fall back to the initial voltage
        BusT plain(Vr0, Vi0);
        plain.allocate();
        plain.initialize();
        plain.evaluateResidual();
        success *= isEqual(plain.Vr(), Vr0);
        success *= isEqual(plain.Vi(), Vi0);

        return success.report(__func__);
      }

      /// Output ports publish the accumulated current injections
      TestOutcome currentOutputs()
      {
        TestStatus success = true;

        auto ir_node = SignalT({.name = "ir", .signal_id = 2});
        auto ii_node = SignalT({.name = "ii", .signal_id = 3});

        BusT bus(1.0, 0.0);
        bus.getPorts().out.template port<SignalOut::ir>().connect(&ir_node);
        bus.getPorts().out.template port<SignalOut::ii>().connect(&ii_node);

        success *= ir_node.assigned();
        success *= !ir_node.linked();

        bus.allocate();
        bus.initialize();
        success *= (bus.verify() == 0);
        success *= ir_node.linked();
        success *= ii_node.linked();
        success *= (ir_node.getVariableIndex() == INVALID_INDEX<IdxT>);
        success *= (ii_node.getVariableIndex() == INVALID_INDEX<IdxT>);

        bus.evaluateResidual();
        success *= isEqual(ir_node.read(), 0.0);
        success *= isEqual(ii_node.read(), 0.0);

        // Two attached components add their injections
        bus.Ir() += -3.7;
        bus.Ii() += 2.4;
        bus.Ir() += 1.3;
        bus.Ii() += -0.8;
        success  *= isEqual(ir_node.read(), -2.4);
        success  *= isEqual(ii_node.read(), 1.6);

        // Re-evaluating resets the sums
        bus.evaluateResidual();
        success *= isEqual(ir_node.read(), 0.0);
        success *= isEqual(ii_node.read(), 0.0);

        return success.report(__func__);
      }

      /// verify() reports connected inputs without a linked source
      TestOutcome verifyUnlinked()
      {
        TestStatus success = true;

        auto vr_node = SignalT({.name = "vr", .signal_id = 0});
        auto vi_node = SignalT({.name = "vi", .signal_id = 1});

        BusT bus(1.0, 0.0);
        bus.getPorts().in.template port<SignalIn::vr>().connect(&vr_node);
        bus.getPorts().in.template port<SignalIn::vi>().connect(&vi_node);
        bus.allocate();
        bus.initialize();

        success *= (bus.verify() == 2);

        ScalarT Vr{0.1};
        IdxT    vr_index{0};
        vr_node.link(&Vr, &vr_index);
        success *= (bus.verify() == 1);

        return success.report(__func__);
      }

      /// Bus voltage carries the dependencies of the input signals
      TestOutcome dependencyTracking()
      {
        TestStatus success = true;

        using VariableT = DependencyTracking::Variable;
        using DtBusT    = PhasorDynamics::BusSignalVoltageIn<VariableT, IdxT>;
        using DtSignalT = PhasorDynamics::SignalNode<VariableT, IdxT>;

        const size_t vr_var_number{10};
        const size_t vi_var_number{11};

        VariableT Vr{0.93};
        VariableT Vi{-0.27};
        Vr.setVariableNumber(vr_var_number);
        Vi.setVariableNumber(vi_var_number);
        IdxT vr_index{5};
        IdxT vi_index{6};

        auto vr_node = DtSignalT({.name = "vr", .signal_id = 0});
        auto vi_node = DtSignalT({.name = "vi", .signal_id = 1});
        vr_node.link(&Vr, &vr_index);
        vi_node.link(&Vi, &vi_index);

        DtBusT bus(VariableT{1.0}, VariableT{0.0});
        bus.getPorts().in.template port<SignalIn::vr>().connect(&vr_node);
        bus.getPorts().in.template port<SignalIn::vi>().connect(&vi_node);
        bus.allocate();
        bus.initialize();
        bus.evaluateResidual();
        bus.evaluateJacobian();

        success *= isEqual(bus.Vr().getValue(), 0.93);
        success *= isEqual(bus.Vi().getValue(), -0.27);

        VariableT::DependencyMap expected_vr{{vr_var_number, 1.0}};
        VariableT::DependencyMap expected_vi{{vi_var_number, 1.0}};
        success *= isEqual(bus.Vr().getDependencies(), expected_vr);
        success *= isEqual(bus.Vi().getDependencies(), expected_vi);

        // A component injection depending on the bus voltage propagates
        bus.Ir() += 2.0 * bus.Vr();
        VariableT::DependencyMap expected_ir{{vr_var_number, 2.0}};
        success *= isEqual(bus.Ir().getDependencies(), expected_ir);

        return success.report(__func__);
      }
    };

  } // namespace Testing
} // namespace GridKit
