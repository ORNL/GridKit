#pragma once

#include <iostream>
#include <map>
#include <stdexcept>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Constants.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusSignalVoltageOut/BusSignalVoltageOut.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusSignalVoltageOut/BusSignalVoltageOutData.hpp>
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
    class BusSignalVoltageOutTests
    {
    public:
      using BusT      = PhasorDynamics::BusSignalVoltageOut<ScalarT, IdxT>;
      using BusTypeT  = typename BusT::BusTypeT;
      using SignalT   = PhasorDynamics::SignalNode<ScalarT, IdxT>;
      using SignalIn  = PhasorDynamics::BusSignalVoltageOutInputs;
      using SignalOut = PhasorDynamics::BusSignalVoltageOutOutputs;

      BusSignalVoltageOutTests()  = default;
      ~BusSignalVoltageOutTests() = default;

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

        const ScalarT Vr{0.93};
        const ScalarT Vi{-0.27};

        PhasorDynamics::BusBase<ScalarT, IdxT>* bus = nullptr;

        bus = new BusT();
        bus->allocate();
        bus->initialize();
        success *= isEqual(bus->Vr(), 0.0);
        success *= isEqual(bus->Vi(), 0.0);
        success *= (bus->size() == 2);
        success *= (bus->BusType() == BusTypeT::SIGNAL_VOLTAGE_OUT);
        delete bus;

        bus = new BusT(Vr, Vi);
        bus->allocate();
        bus->initialize();
        success *= isEqual(bus->Vr(), Vr);
        success *= isEqual(bus->Vi(), Vi);
        delete bus;

        bus = nullptr;

        return success.report(__func__);
      }

      /// Signal outlets publish the bus voltage and its variable indices
      TestOutcome voltageOutputs()
      {
        TestStatus success = true;

        const ScalarT Vr{0.93};
        const ScalarT Vi{-0.27};
        const IdxT    vr_index{7};
        const IdxT    vi_index{8};

        auto vr_node = SignalT({.name = "vr", .signal_id = 0});
        auto vi_node = SignalT({.name = "vi", .signal_id = 1});

        // Mandatory current inlets
        ScalarT Ir{-3.7};
        ScalarT Ii{2.4};
        IdxT    ir_index{5};
        IdxT    ii_index{6};
        auto    ir_node = SignalT({.name = "ir", .signal_id = 2});
        auto    ii_node = SignalT({.name = "ii", .signal_id = 3});
        ir_node.link(&Ir, &ir_index);
        ii_node.link(&Ii, &ii_index);

        BusT bus(Vr, Vi);
        bus.getPorts().out.template port<SignalOut::vr>().connect(&vr_node);
        bus.getPorts().out.template port<SignalOut::vi>().connect(&vi_node);
        bus.getPorts().in.template port<SignalIn::ir>().connect(&ir_node);
        bus.getPorts().in.template port<SignalIn::ii>().connect(&ii_node);

        success *= vr_node.assigned();
        success *= vi_node.assigned();
        success *= !vr_node.linked();

        bus.allocate();
        bus.setVariableIndex(0, vr_index);
        bus.setVariableIndex(1, vi_index);
        bus.initialize();

        success *= (bus.verify() == 0);
        success *= vr_node.linked();
        success *= vi_node.linked();
        success *= isEqual(vr_node.read(), Vr);
        success *= isEqual(vi_node.read(), Vi);
        success *= (vr_node.getVariableIndex() == vr_index);
        success *= (vi_node.getVariableIndex() == vi_index);

        // Signal follows the live bus voltage
        bus.Vr()  = 1.17;
        bus.Vi()  = 0.41;
        success  *= isEqual(vr_node.read(), 1.17);
        success  *= isEqual(vi_node.read(), 0.41);

        return success.report(__func__);
      }

      /// Signal inlets add current injections to the residual
      TestOutcome residual()
      {
        TestStatus success = true;

        const ScalarT Vr{0.93};
        const ScalarT Vi{-0.27};
        ScalarT       Ir{-3.7}; ///< Current injection on signal ir
        ScalarT       Ii{2.4};  ///< Current injection on signal ii
        IdxT          ir_index{5};
        IdxT          ii_index{6};

        // A bus without connected current inlets is rejected by verify()
        {
          BusT bus(Vr, Vi);
          bus.allocate();
          bus.initialize();
          const auto previous_verbosity = Log::verbosity();
          Log::setVerbosity(Log::Verbosity::NONE); // expected errors
          success *= verifyThrows(bus);
          Log::setVerbosity(previous_verbosity);
        }

        auto ir_node = SignalT({.name = "ir", .signal_id = 2});
        auto ii_node = SignalT({.name = "ii", .signal_id = 3});
        ir_node.link(&Ir, &ir_index);
        ii_node.link(&Ii, &ii_index);

        BusT bus(Vr, Vi);
        bus.getPorts().in.template port<SignalIn::ir>().connect(&ir_node);
        bus.getPorts().in.template port<SignalIn::ii>().connect(&ii_node);
        bus.allocate();
        bus.initialize();
        success *= (bus.verify() == 0);

        bus.evaluateResidual();
        success *= isEqual(bus.Ir(), Ir);
        success *= isEqual(bus.Ii(), Ii);
        success *= isEqual(bus.getResidual().getData()[0], Ir);
        success *= isEqual(bus.getResidual().getData()[1], Ii);

        // A device attached to the bus adds its current after the bus residual
        bus.Ir() += 1.3;
        bus.Ii() += -0.8;
        success  *= isEqual(bus.Ir(), Ir + 1.3);
        success  *= isEqual(bus.Ii(), Ii - 0.8);

        // Re-evaluating resets and re-reads the signals
        Ir = 0.55;
        Ii = -1.25;
        bus.evaluateResidual();
        success *= isEqual(bus.Ir(), Ir);
        success *= isEqual(bus.Ii(), Ii);

        return success.report(__func__);
      }

      /// verify() throws for current inlets that are unconnected or unlinked
      TestOutcome verifyUnlinked()
      {
        TestStatus success = true;

        // This test triggers error messages on purpose; silence them.
        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::NONE);

        auto ir_node = SignalT({.name = "ir", .signal_id = 2});
        auto ii_node = SignalT({.name = "ii", .signal_id = 3});

        BusT bus(1.0, 0.0);
        bus.getPorts().in.template port<SignalIn::ir>().connect(&ir_node);
        bus.getPorts().in.template port<SignalIn::ii>().connect(&ii_node);
        bus.allocate();
        bus.initialize();

        success *= verifyThrows(bus);

        ScalarT Ir{0.1};
        IdxT    ir_index{0};
        ir_node.link(&Ir, &ir_index);
        success *= verifyThrows(bus);

        ScalarT Ii{0.2};
        IdxT    ii_index{1};
        ii_node.link(&Ii, &ii_index);
        success *= (bus.verify() == 0);

        Log::setVerbosity(previous_verbosity);

        return success.report(__func__);
      }

      /// Residual dependencies on signal variables via dependency tracking
      TestOutcome dependencyTracking()
      {
        TestStatus success = true;

        using VariableT = DependencyTracking::Variable;
        using DtBusT    = PhasorDynamics::BusSignalVoltageOut<VariableT, IdxT>;
        using DtSignalT = PhasorDynamics::SignalNode<VariableT, IdxT>;

        const size_t ir_var_number{10};
        const size_t ii_var_number{11};

        VariableT Ir{-3.7};
        VariableT Ii{2.4};
        Ir.setVariableNumber(ir_var_number);
        Ii.setVariableNumber(ii_var_number);
        IdxT ir_index{5};
        IdxT ii_index{6};

        auto ir_node = DtSignalT({.name = "ir", .signal_id = 2});
        auto ii_node = DtSignalT({.name = "ii", .signal_id = 3});
        ir_node.link(&Ir, &ir_index);
        ii_node.link(&Ii, &ii_index);

        DtBusT bus(VariableT{0.93}, VariableT{-0.27});
        bus.getPorts().in.template port<SignalIn::ir>().connect(&ir_node);
        bus.getPorts().in.template port<SignalIn::ii>().connect(&ii_node);
        bus.allocate();
        bus.setVariableIndex(0, 3);
        bus.setVariableIndex(1, 4);
        bus.initialize();
        bus.evaluateResidual();
        bus.evaluateJacobian();

        const auto* f = bus.getResidual().getData();

        success *= isEqual(f[0].getValue(), -3.7);
        success *= isEqual(f[1].getValue(), 2.4);

        VariableT::DependencyMap expected_f0{{ir_var_number, 1.0}};
        VariableT::DependencyMap expected_f1{{ii_var_number, 1.0}};
        success *= isEqual(f[0].getDependencies(), expected_f0);
        success *= isEqual(f[1].getDependencies(), expected_f1);

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /// Sparse Jacobian has unit entries in the signal columns
      TestOutcome jacobian()
      {
        TestStatus success = true;

        ScalarT    Ir{-3.7};
        ScalarT    Ii{2.4};
        IdxT       ir_index{5};
        IdxT       ii_index{6};
        const IdxT var_offset{3};

        auto ir_node = SignalT({.name = "ir", .signal_id = 2});
        auto ii_node = SignalT({.name = "ii", .signal_id = 3});
        ir_node.link(&Ir, &ir_index);
        ii_node.link(&Ii, &ii_index);

        BusT bus(0.93, -0.27);
        bus.getPorts().in.template port<SignalIn::ir>().connect(&ir_node);
        bus.getPorts().in.template port<SignalIn::ii>().connect(&ii_node);
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
        success *= (jac->getNnz() == 6);

        // Accumulate COO triplets into per-row maps (duplicates are summed)
        std::map<size_t, std::map<size_t, double>> rows;
        const auto*                                row_data = jac->getRowData();
        const auto*                                col_data = jac->getColData();
        const auto*                                values   = jac->getValues();
        for (IdxT k = 0; k < jac->getNnz(); ++k)
        {
          rows[static_cast<size_t>(row_data[k])][static_cast<size_t>(col_data[k])] += values[k];
        }

        std::map<size_t, double> expected_row0{{3, 0.0}, {4, 0.0}, {5, 1.0}};
        std::map<size_t, double> expected_row1{{3, 0.0}, {4, 0.0}, {6, 1.0}};
        success *= (rows.size() == 2);
        success *= isEqual(rows[3], expected_row0);
        success *= isEqual(rows[4], expected_row1);

        return success.report(__func__);
      }
#endif
    };

  } // namespace Testing
} // namespace GridKit
