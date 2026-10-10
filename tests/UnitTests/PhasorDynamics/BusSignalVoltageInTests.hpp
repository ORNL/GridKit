#pragma once

#include <iostream>
#include <stdexcept>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Constants.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusSignalVoltageIn/BusSignalVoltageIn.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeData.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace Testing
  {
    using Log = ::GridKit::Utilities::Logger;

    template <class ScalarT, typename IdxT>
    class BusSignalVoltageInTests
    {
    public:
      using BusT      = PhasorDynamics::BusSignalVoltageIn<ScalarT, IdxT>;
      using BusTypeT  = typename BusT::BusTypeT;
      using SignalT   = PhasorDynamics::SignalNode<ScalarT, IdxT>;
      using SignalIn  = PhasorDynamics::BusSignalInputs;
      using SignalOut = PhasorDynamics::BusSignalOutputs;

      BusSignalVoltageInTests()  = default;
      ~BusSignalVoltageInTests() = default;

      /// Invalid configurations are reported without throwing.
      template <typename BusLike>
      static bool verifyFails(const BusLike& bus)
      {
        return !bus.verify().passed();
      }

      /// Constructor, allocation, and initialization checks
      TestOutcome constructor()
      {
        TestStatus success = true;

        // Keep expected invalid-configuration diagnostics quiet.
        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::NONE);

        const ScalarT Vr{0.93};
        const ScalarT Vi{-0.27};

        PhasorDynamics::BusBase<ScalarT, IdxT>* bus = nullptr;

        bus = new BusT();
        bus->allocate();
        bus->initialize();
        success *= (bus->size() == 0);
        success *= (bus->BusType() == BusTypeT::SIGNAL_VOLTAGE_IN);
        success *= isEqual(bus->Ir(), 0.0);
        success *= isEqual(bus->Ii(), 0.0);
        // Voltage inlets are mandatory: an unconnected bus fails verification
        success *= verifyFails(*bus);
        delete bus;

        // Initial voltage arguments are accepted for interface uniformity but not used
        bus = new BusT(Vr, Vi);
        bus->allocate();
        bus->initialize();
        success *= verifyFails(*bus);
        delete bus;

        bus = nullptr;

        Log::setVerbosity(previous_verbosity);

        return success.report(__func__);
      }

      /// Signal inlets set the bus voltage
      TestOutcome voltageInputs()
      {
        TestStatus success = true;

        // Keep expected invalid-configuration diagnostics quiet.
        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::NONE);

        ScalarT Vr{0.93};  ///< Voltage on signal vr
        ScalarT Vi{-0.27}; ///< Voltage on signal vi
        IdxT    vr_index{7};
        IdxT    vi_index{8};

        auto vr_node = SignalT({.name = "vr", .signal_id = 0});
        auto vi_node = SignalT({.name = "vi", .signal_id = 1});
        vr_node.link(&Vr, &vr_index);
        vi_node.link(&Vi, &vi_index);

        BusT bus;
        bus.getPorts().in.template port<SignalIn::vr>().connect(&vr_node);
        bus.getPorts().in.template port<SignalIn::vi>().connect(&vi_node);
        bus.allocate();
        bus.initialize();
        success *= bus.verify().passed();

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

        // Reading an unconnected voltage inlet is an error, never a default value
        BusT plain;
        plain.allocate();
        plain.initialize();
        success    *= verifyFails(plain);
        bool threw  = false;
        try
        {
          [[maybe_unused]] const auto& v = plain.Vr();
        }
        catch (const std::runtime_error&)
        {
          threw = true;
        }
        success *= threw;

        Log::setVerbosity(previous_verbosity);

        return success.report(__func__);
      }

      /// Signal outlets publish the accumulated current injections
      TestOutcome currentOutputs()
      {
        TestStatus success = true;

        // Mandatory voltage inlets
        ScalarT Vr{0.93};
        ScalarT Vi{-0.27};
        IdxT    vr_index{7};
        IdxT    vi_index{8};
        auto    vr_node = SignalT({.name = "vr", .signal_id = 0});
        auto    vi_node = SignalT({.name = "vi", .signal_id = 1});
        vr_node.link(&Vr, &vr_index);
        vi_node.link(&Vi, &vi_index);

        auto ir_node = SignalT({.name = "ir", .signal_id = 2});
        auto ii_node = SignalT({.name = "ii", .signal_id = 3});

        BusT bus;
        bus.getPorts().in.template port<SignalIn::vr>().connect(&vr_node);
        bus.getPorts().in.template port<SignalIn::vi>().connect(&vi_node);
        bus.getPorts().out.template port<SignalOut::ir>().connect(&ir_node);
        bus.getPorts().out.template port<SignalOut::ii>().connect(&ii_node);

        success *= ir_node.assigned();
        success *= !ir_node.linked();

        bus.allocate();
        bus.initialize();
        success *= bus.verify().passed();
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

      /// verify() reports errors for voltage inlets that are unconnected or unlinked
      TestOutcome verifyUnlinked()
      {
        TestStatus success = true;

        // Keep expected invalid-configuration diagnostics quiet.
        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::NONE);

        auto vr_node = SignalT({.name = "vr", .signal_id = 0});
        auto vi_node = SignalT({.name = "vi", .signal_id = 1});

        BusT bus;
        bus.getPorts().in.template port<SignalIn::vr>().connect(&vr_node);
        bus.allocate();
        bus.initialize();

        // vr connected but unlinked, vi not connected
        success *= verifyFails(bus);

        ScalarT Vr{0.1};
        IdxT    vr_index{0};
        vr_node.link(&Vr, &vr_index);
        success *= verifyFails(bus);

        bus.getPorts().in.template port<SignalIn::vi>().connect(&vi_node);
        success *= verifyFails(bus);

        ScalarT Vi{0.2};
        IdxT    vi_index{1};
        vi_node.link(&Vi, &vi_index);
        success *= bus.verify().passed();

        Log::setVerbosity(previous_verbosity);

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

        DtBusT bus;
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
