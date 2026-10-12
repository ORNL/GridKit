#pragma once

#include <iostream>
#include <map>
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

      /// True if verify() throws, as it must for a misconnected bus
      template <typename BusLike>
      static bool verifyThrows(const BusLike& bus)
      {
        try
        {
          bus.verify();
        }
        catch (const std::runtime_error&)
        {
          return true;
        }
        return false;
      }

      /// Constructor, allocation, and initialization checks
      TestOutcome constructor()
      {
        TestStatus success = true;

        // This test triggers error messages on purpose; silence them.
        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::NONE);

        const ScalarT Vr{0.93};
        const ScalarT Vi{-0.27};

        PhasorDynamics::BusBase<ScalarT, IdxT>* bus = nullptr;

        bus = new BusT();
        bus->allocate();
        bus->initialize();
        success *= (bus->size() == 2);
        success *= (bus->BusType() == BusTypeT::SIGNAL_VOLTAGE_IN);
        success *= isEqual(bus->Ir(), 0.0);
        success *= isEqual(bus->Ii(), 0.0);
        // Voltage inlets are mandatory: an unconnected bus fails verification
        success *= verifyThrows(*bus);
        delete bus;

        // Initial voltage arguments are accepted for interface uniformity but not used
        bus = new BusT(Vr, Vi);
        bus->allocate();
        bus->initialize();
        success *= verifyThrows(*bus);
        delete bus;

        bus = nullptr;

        Log::setVerbosity(previous_verbosity);

        return success.report(__func__);
      }

      /// Signal inlets set the bus voltage
      TestOutcome voltageInputs()
      {
        TestStatus success = true;

        // This test triggers error messages on purpose; silence them.
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

        // Terminal nodes the bus publishes its voltage on for attached components
        SignalT terminal_vr;
        SignalT terminal_vi;

        BusT bus;
        bus.getPorts().in.template port<SignalIn::vr>().connect(&vr_node);
        bus.getPorts().in.template port<SignalIn::vi>().connect(&vi_node);
        bus.getPorts().out.template port<SignalOut::vr>().connect(&terminal_vr);
        bus.getPorts().out.template port<SignalOut::vi>().connect(&terminal_vi);
        bus.allocate();
        bus.initialize();
        success *= (bus.verify() == 0);

        // Terminal nodes carry the inlet voltage and column, and the bus row
        success *= terminal_vr.linked();
        success *= terminal_vi.linked();
        success *= (&terminal_vr.read() == &Vr);
        success *= (&terminal_vi.read() == &Vi);
        success *= (terminal_vr.getVariableIndex() == vr_index);
        success *= (terminal_vi.getVariableIndex() == vi_index);
        success *= (terminal_vr.getResidualIndex() == bus.getResidualIndex(0));
        success *= (terminal_vi.getResidualIndex() == bus.getResidualIndex(1));

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
        success    *= verifyThrows(plain);
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

      /// Signal outlets publish the boundary currents, which are the bus variables
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
        success *= (bus.verify() == 0);
        success *= ir_node.linked();
        success *= ii_node.linked();
        success *= (ir_node.getVariableIndex() == bus.getVariableIndex(0));
        success *= (ii_node.getVariableIndex() == bus.getVariableIndex(1));
        success *= isEqual(ir_node.read(), 0.0);
        success *= isEqual(ii_node.read(), 0.0);

        // The outlets publish the boundary current variables
        bus.Ir()  = -2.4;
        bus.Ii()  = 1.6;
        success  *= isEqual(ir_node.read(), -2.4);
        success  *= isEqual(ii_node.read(), 1.6);

        // Each residual row is the injection sum minus the boundary current;
        // with nothing attached the sum is zero.
        bus.evaluateResidual();
        const auto* f  = bus.getResidual().getData();
        success       *= isEqual(f[0], 2.4);
        success       *= isEqual(f[1], -1.6);

        return success.report(__func__);
      }

      /// verify() throws for voltage inlets that are unconnected or unlinked
      TestOutcome verifyUnlinked()
      {
        TestStatus success = true;

        // This test triggers error messages on purpose; silence them.
        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::NONE);

        auto vr_node = SignalT({.name = "vr", .signal_id = 0});
        auto vi_node = SignalT({.name = "vi", .signal_id = 1});

        BusT bus;
        bus.getPorts().in.template port<SignalIn::vr>().connect(&vr_node);
        bus.allocate();
        bus.initialize();

        // vr connected but unlinked, vi not connected
        success *= verifyThrows(bus);

        ScalarT Vr{0.1};
        IdxT    vr_index{0};
        vr_node.link(&Vr, &vr_index);
        success *= verifyThrows(bus);

        bus.getPorts().in.template port<SignalIn::vi>().connect(&vi_node);
        success *= verifyThrows(bus);

        ScalarT Vi{0.2};
        IdxT    vi_index{1};
        vi_node.link(&Vi, &vi_index);
        success *= (bus.verify() == 0);

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

        // Each residual row depends only on its own boundary current
        const auto*              f = bus.getResidual().getData();
        VariableT::DependencyMap expected_f0{{2 * static_cast<size_t>(bus.getVariableIndex(0)), -1.0}};
        VariableT::DependencyMap expected_f1{{2 * static_cast<size_t>(bus.getVariableIndex(1)), -1.0}};
        success *= isEqual(f[0].getDependencies(), expected_f0);
        success *= isEqual(f[1].getDependencies(), expected_f1);

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /// The Jacobian is minus identity on the boundary currents
      TestOutcome jacobian()
      {
        TestStatus success = true;

        ScalarT    Vr{0.93};
        ScalarT    Vi{-0.27};
        IdxT       vr_index{7};
        IdxT       vi_index{8};
        const IdxT var_offset{3};

        auto vr_node = SignalT({.name = "vr", .signal_id = 0});
        auto vi_node = SignalT({.name = "vi", .signal_id = 1});
        vr_node.link(&Vr, &vr_index);
        vi_node.link(&Vi, &vi_index);

        BusT bus;
        bus.getPorts().in.template port<SignalIn::vr>().connect(&vr_node);
        bus.getPorts().in.template port<SignalIn::vi>().connect(&vi_node);
        bus.allocate();
        for (IdxT i = 0; i < bus.size(); ++i)
        {
          bus.setVariableIndex(i, i + var_offset);
          bus.setResidualIndex(i, i + var_offset);
        }
        bus.initialize();
        bus.evaluateResidual();
        bus.evaluateJacobian();

        auto* jac  = bus.getCooJacobian();
        success   *= (jac != nullptr);
        if (jac == nullptr)
        {
          return success.report(__func__);
        }

        // Accumulate COO triplets into per-row maps (duplicates are summed)
        std::map<size_t, std::map<size_t, double>> rows;
        const auto*                                row_data = jac->getRowData();
        const auto*                                col_data = jac->getColData();
        const auto*                                values   = jac->getValues();
        for (IdxT k = 0; k < jac->getNnz(); ++k)
        {
          rows[static_cast<size_t>(row_data[k])][static_cast<size_t>(col_data[k])] += values[k];
        }

        std::map<size_t, double> expected_row0{{3, -1.0}};
        std::map<size_t, double> expected_row1{{4, -1.0}};
        success *= (rows.size() == 2);
        success *= isEqual(rows[3], expected_row0);
        success *= isEqual(rows[4], expected_row1);

        return success.report(__func__);
      }
#endif
    };

  } // namespace Testing
} // namespace GridKit
