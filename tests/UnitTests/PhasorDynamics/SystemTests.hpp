#pragma once

#include <cstddef>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Branch/Branch.hpp>
#include <GridKit/Model/PhasorDynamics/Branch/BranchData.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusSignalVoltageIn/BusSignalVoltageIn.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusSignalVoltageOut/BusSignalVoltageOut.hpp>
#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/Load/LoadZ/LoadZ.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    using GridKit::PhasorDynamics::BranchBuses;
    using GridKit::PhasorDynamics::BranchParameters;

    using Log = ::GridKit::Utilities::Logger;

    template <class ScalarT, typename IdxT>
    class SystemTests
    {
    private:
      using ComponentT = PhasorDynamics::Component<ScalarT, IdxT>;
      using RealT      = typename ComponentT::RealT;

      class InitializationFailureComponent final : public ComponentT
      {
      public:
        InitializationFailureComponent()
        {
          this->size_ = static_cast<IdxT>(1);
        }

        int setGridKitComponentID(IdxT component_id) override final
        {
          this->gridkit_component_id_ = component_id;
          return 0;
        }

        int allocate() override final
        {
          if (!this->allocated_)
          {
            this->allocateVectors(this->size_);
          }

          const auto size = static_cast<std::size_t>(this->size_);
          this->tag_.assign(size, false);
          this->variable_indices_.resize(size);
          this->residual_indices_.resize(size);
          this->allocated_ = true;
          return 0;
        }

        int verify() const override final
        {
          return 0;
        }

        int initialize() override final
        {
          return 1;
        }

        int tagDifferentiable() override final
        {
          return 0;
        }

        int setAbsoluteTolerance(RealT) override final
        {
          return 0;
        }

        int evaluateResidual() override final
        {
          return 0;
        }

        int evaluateJacobian() override final
        {
          return this->constructCoo();
        }
      };

      /// Zero-size component that exposes the system base it was given
      class BaseProbeComponent final : public ComponentT
      {
      public:
        int setGridKitComponentID(IdxT component_id) override final
        {
          this->gridkit_component_id_ = component_id;
          return 0;
        }

        int allocate() override final
        {
          this->allocated_ = true;
          return 0;
        }

        int verify() const override final
        {
          return 0;
        }

        int initialize() override final
        {
          return 0;
        }

        int tagDifferentiable() override final
        {
          return 0;
        }

        int setAbsoluteTolerance(RealT) override final
        {
          return 0;
        }

        int evaluateResidual() override final
        {
          return 0;
        }

        int evaluateJacobian() override final
        {
          return this->constructCoo();
        }

        RealT vaSystemBase() const
        {
          return this->va_system_base_;
        }
      };

    public:
      SystemTests()  = default;
      ~SystemTests() = default;

      /// Constructor, allocation, and initialization checks
      TestOutcome constructor()
      {
        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT>* system = nullptr;

        // Create an empty system
        system = new PhasorDynamics::SystemModel<ScalarT, IdxT>();

        if (system == nullptr)
        {
          std::cout << "Default constructor failed!\n";
          success = false;
          return success.report(__func__);
        }

        delete system;
        system = nullptr;

        PhasorDynamics::SystemModelData<ScalarT, IdxT> data;

        // Set bus data
        data.bus.resize(2);

        // Bus 0
        data.bus[0].bus_id   = 0;
        data.bus[0].bus_type = PhasorDynamics::BusData<ScalarT, IdxT>::BusType::SLACK;
        data.bus[0].Vr0      = 10.0;
        data.bus[0].Vi0      = 20.0;

        // Bus 1
        data.bus[1].bus_id   = 1;
        data.bus[1].bus_type = PhasorDynamics::BusData<ScalarT, IdxT>::BusType::DEFAULT;
        data.bus[1].Vr0      = 30.0;
        data.bus[1].Vi0      = 40.0;

        // Set branch data
        data.branch.resize(1);

        // Branch 0-1
        data.branch[0].buses[BranchBuses::bus1]        = data.bus[0].bus_id;
        data.branch[0].buses[BranchBuses::bus2]        = data.bus[1].bus_id;
        data.branch[0].parameters[BranchParameters::R] = 2.0;
        data.branch[0].parameters[BranchParameters::X] = 4.0;
        data.branch[0].parameters[BranchParameters::G] = 0.2;
        data.branch[0].parameters[BranchParameters::B] = 1.2;

        // Create an empty system model
        system = new PhasorDynamics::SystemModel<ScalarT, IdxT>(data);
        system->allocate();
        system->initialize();
        system->evaluateResidual();

        // Answer keys
        const ScalarT Ir0{17.0};  ///< Solution: real current entering bus-0
        const ScalarT Ii0{-10.0}; ///< Solution: imaginary current entering bus-0
        const ScalarT Ir1{15.0};  ///< Solution: real current entering bus-1
        const ScalarT Ii1{-20.0}; ///< Solution: imaginary current entering bus-1

        // Components are stored in insertion order and the composer adds the
        // two buses before the branch.
        auto* branch = system->getComponent(2);
        auto* bus1   = system->getBus(1);

        success *= isEqual(branch->getExternalResidual()[0], Ir0);
        success *= isEqual(branch->getExternalResidual()[1], Ii0);
        success *= isEqual(bus1->Ir(), Ir1);
        success *= isEqual(bus1->Ii(), Ii1);

        delete system;
        system = nullptr;

        return success.report(__func__);
      }

      TestOutcome composer()
      {
        TestStatus success = true;

        RealT R{2.0}; ///< Branch series resistance
        RealT X{4.0}; ///< Branch series reactance
        RealT G{0.2}; ///< Branch shunt conductance
        RealT B{1.2}; ///< Branch shunt charging

        ScalarT Vr1{10.0}; ///< Bus-1 real voltage
        ScalarT Vi1{20.0}; ///< Bus-1 imaginary voltage
        ScalarT Vr2{30.0}; ///< Bus-2 real voltage
        ScalarT Vi2{40.0}; ///< Bus-2 imaginary voltage

        const ScalarT Ir1{17.0};  ///< Solution: real current entering bus-1
        const ScalarT Ii1{-10.0}; ///< Solution: imaginary current entering bus-1
        const ScalarT Ir2{15.0};  ///< Solution: real current entering bus-2
        const ScalarT Ii2{-20.0}; ///< Solution: imaginary current entering bus-2

        // Create an empty system model
        PhasorDynamics::SystemModel<ScalarT, IdxT> system;

        // Add a bus
        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus1(Vr1, Vi1);
        system.addBus(&bus1);

        // Add a bus
        PhasorDynamics::Bus<ScalarT, IdxT> bus2(Vr2, Vi2);
        system.addBus(&bus2);

        PhasorDynamics::SignalNode<ScalarT, IdxT> bus1_vr;
        PhasorDynamics::SignalNode<ScalarT, IdxT> bus1_vi;
        PhasorDynamics::SignalNode<ScalarT, IdxT> bus2_vr;
        PhasorDynamics::SignalNode<ScalarT, IdxT> bus2_vi;
        bus1.getPorts().out[PhasorDynamics::BusSignalOutputs::vr].connect(&bus1_vr);
        bus1.getPorts().out[PhasorDynamics::BusSignalOutputs::vi].connect(&bus1_vi);
        bus2.getPorts().out[PhasorDynamics::BusSignalOutputs::vr].connect(&bus2_vr);
        bus2.getPorts().out[PhasorDynamics::BusSignalOutputs::vi].connect(&bus2_vi);

        PhasorDynamics::Branch<ScalarT, IdxT> branch(R, X, G, B);
        branch.getPorts().terminal_vr[PhasorDynamics::BranchBuses::bus1].connect(&bus1_vr);
        branch.getPorts().terminal_vi[PhasorDynamics::BranchBuses::bus1].connect(&bus1_vi);
        branch.getPorts().terminal_vr[PhasorDynamics::BranchBuses::bus2].connect(&bus2_vr);
        branch.getPorts().terminal_vi[PhasorDynamics::BranchBuses::bus2].connect(&bus2_vi);
        system.addComponent(&branch);

        system.allocate();
        system.initialize();
        system.evaluateResidual();

        success *= isEqual(branch.getExternalResidual()[0], Ir1);
        success *= isEqual(branch.getExternalResidual()[1], Ii1);
        success *= isEqual(bus2.Ir(), Ir2);
        success *= isEqual(bus2.Ii(), Ii2);

        system.evaluateResidual();

        success *= isEqual(branch.getExternalResidual()[0], Ir1);
        success *= isEqual(branch.getExternalResidual()[1], Ii1);
        success *= isEqual(bus2.Ir(), Ir2);
        success *= isEqual(bus2.Ii(), Ii2);

        return success.report(__func__);
      }

      /// A load on a bus whose voltage comes from signals: the system adds the
      /// load current into the bus rows, and the `ir`/`ii` outlets publish the
      /// boundary currents as system variables.
      TestOutcome signalVoltageInAssembly()
      {
        TestStatus success = true;

        using SignalIn  = PhasorDynamics::BusSignalInputs;
        using SignalOut = PhasorDynamics::BusSignalOutputs;

        // Voltage arrives from outside the system, so it has no column.
        ScalarT Vr{1.06};
        ScalarT Vi{-0.01};
        IdxT    vr_index{INVALID_INDEX<IdxT>};
        IdxT    vi_index{INVALID_INDEX<IdxT>};

        PhasorDynamics::SignalNode<ScalarT, IdxT> vr_inlet;
        PhasorDynamics::SignalNode<ScalarT, IdxT> vi_inlet;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ir_outlet;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ii_outlet;
        PhasorDynamics::SignalNode<ScalarT, IdxT> bus_vr;
        PhasorDynamics::SignalNode<ScalarT, IdxT> bus_vi;
        vr_inlet.link(&Vr, &vr_index);
        vi_inlet.link(&Vi, &vi_index);

        PhasorDynamics::SystemModel<ScalarT, IdxT>        system;
        PhasorDynamics::BusSignalVoltageIn<ScalarT, IdxT> bus;
        PhasorDynamics::LoadZ<ScalarT, IdxT>              load(0.1, 0.0);

        bus.getPorts().in.template port<SignalIn::vr>().connect(&vr_inlet);
        bus.getPorts().in.template port<SignalIn::vi>().connect(&vi_inlet);
        bus.getPorts().out.template port<SignalOut::ir>().connect(&ir_outlet);
        bus.getPorts().out.template port<SignalOut::ii>().connect(&ii_outlet);
        bus.getPorts().out[PhasorDynamics::BusSignalOutputs::vr].connect(&bus_vr);
        bus.getPorts().out[PhasorDynamics::BusSignalOutputs::vi].connect(&bus_vi);
        load.getPorts().terminal_vr[PhasorDynamics::LoadZBuses::bus].connect(&bus_vr);
        load.getPorts().terminal_vi[PhasorDynamics::LoadZBuses::bus].connect(&bus_vi);

        system.addBus(&bus);
        system.addComponent(&load);

        success *= system.allocate() == 0;
        success *= system.initialize() == 0;
        success *= system.evaluateResidual() == 0;

        success *= ir_outlet.getVariableIndex() == bus.getVariableIndex(0);
        success *= ii_outlet.getVariableIndex() == bus.getVariableIndex(1);

        // With zero boundary current, each bus row holds the load injection.
        const auto* f       = system.getResidual().getData();
        const auto  ir_row  = static_cast<size_t>(bus.getResidualIndex(0));
        const auto  ii_row  = static_cast<size_t>(bus.getResidualIndex(1));
        const auto  Ir      = load.getExternalResidual()[0];
        const auto  Ii      = load.getExternalResidual()[1];
        success            *= !isEqual(Ir, 0.0);
        success            *= isEqual(f[ir_row], Ir);
        success            *= isEqual(f[ii_row], Ii);

        // The consistent boundary current zeroes the rows and is published.
        bus.Ir()  = Ir;
        bus.Ii()  = Ii;
        success  *= system.evaluateResidual() == 0;
        success  *= isEqual(f[ir_row], 0.0);
        success  *= isEqual(f[ii_row], 0.0);
        success  *= isEqual(ir_outlet.read(), Ir);
        success  *= isEqual(ii_outlet.read(), Ii);

        return success.report(__func__);
      }

      /// Global index writes on the system route to the owning component
      TestOutcome indexRouting()
      {
        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT> system;

        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus1(10.0, 20.0);
        system.addBus(&bus1);

        PhasorDynamics::Bus<ScalarT, IdxT> bus2(30.0, 40.0);
        system.addBus(&bus2);

        PhasorDynamics::SignalNode<ScalarT, IdxT> bus1_vr;
        PhasorDynamics::SignalNode<ScalarT, IdxT> bus1_vi;
        PhasorDynamics::SignalNode<ScalarT, IdxT> bus2_vr;
        PhasorDynamics::SignalNode<ScalarT, IdxT> bus2_vi;
        bus1.getPorts().out[PhasorDynamics::BusSignalOutputs::vr].connect(&bus1_vr);
        bus1.getPorts().out[PhasorDynamics::BusSignalOutputs::vi].connect(&bus1_vi);
        bus2.getPorts().out[PhasorDynamics::BusSignalOutputs::vr].connect(&bus2_vr);
        bus2.getPorts().out[PhasorDynamics::BusSignalOutputs::vi].connect(&bus2_vi);

        PhasorDynamics::Branch<ScalarT, IdxT> branch(2.0, 4.0, 0.2, 1.2);
        branch.getPorts().terminal_vr[PhasorDynamics::BranchBuses::bus1].connect(&bus1_vr);
        branch.getPorts().terminal_vi[PhasorDynamics::BranchBuses::bus1].connect(&bus1_vi);
        branch.getPorts().terminal_vr[PhasorDynamics::BranchBuses::bus2].connect(&bus2_vr);
        branch.getPorts().terminal_vi[PhasorDynamics::BranchBuses::bus2].connect(&bus2_vi);
        system.addComponent(&branch);

        system.allocate();

        // Bus-2 owns both system variables and allocate() default-maps them
        // to the local index.
        success *= (bus2.getVariableIndex(0) == 0);
        success *= (bus2_vr.getVariableIndex() == 0);

        // A parent model reassigns global indices through the system setters.
        for (IdxT j = 0; j < system.size(); ++j)
        {
          system.setVariableIndex(j, j + 100);
          system.setResidualIndex(j, j + 200);
        }

        // Routed writes land in the owning component's index map ...
        success *= (bus2.getVariableIndex(0) == 100);
        success *= (bus2.getVariableIndex(1) == 101);
        success *= (bus2.getResidualIndex(0) == 200);
        success *= (bus2.getResidualIndex(1) == 201);

        // ... where the published signal nodes point ...
        success *= (bus2_vr.getVariableIndex() == 100);
        success *= (bus2_vi.getVariableIndex() == 101);
        success *= (bus2_vr.getResidualIndex() == 200);
        success *= (bus2_vi.getResidualIndex() == 201);

        // ... and consumers gather them through their attached signals.
        success *= (branch.getPorts().terminal_vr[PhasorDynamics::BranchBuses::bus2].signalVariableIndex() == 100);
        success *= (branch.getPorts().terminal_vi[PhasorDynamics::BranchBuses::bus2].signalVariableIndex() == 101);
        success *= (branch.getPorts().terminal_vr[PhasorDynamics::BranchBuses::bus2].signalResidualIndex() == 200);
        success *= (branch.getPorts().terminal_vi[PhasorDynamics::BranchBuses::bus2].signalResidualIndex() == 201);

        // The system's own index map mirrors the assignment.
        success *= (system.getVariableIndex(0) == 100);
        success *= (system.getResidualIndex(1) == 201);

        // Constant slack voltages carry no solver indices before or after.
        success *= (bus1_vr.getVariableIndex() == INVALID_INDEX<IdxT>);
        success *= (branch.getPorts().terminal_vr[PhasorDynamics::BranchBuses::bus1].signalVariableIndex() == INVALID_INDEX<IdxT>);

        return success.report(__func__);
      }

      /// A system model binds into a parent system like any other component
      TestOutcome nestedSystem()
      {
        TestStatus success = true;

        RealT R{2.0}; ///< Branch series resistance
        RealT X{4.0}; ///< Branch series reactance
        RealT G{0.2}; ///< Branch shunt conductance
        RealT B{1.2}; ///< Branch shunt charging

        const ScalarT Ir1{17.0};  ///< Solution: real current entering bus-1
        const ScalarT Ii1{-10.0}; ///< Solution: imaginary current entering bus-1
        const ScalarT Ir2{15.0};  ///< Solution: real current entering bus-2
        const ScalarT Ii2{-20.0}; ///< Solution: imaginary current entering bus-2

        // The root system owns the slack bus; a nested system owns the PQ
        // bus and the branch.
        PhasorDynamics::SystemModel<ScalarT, IdxT> root;
        PhasorDynamics::SystemModel<ScalarT, IdxT> nested;

        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus1(10.0, 20.0);
        root.addBus(&bus1);

        PhasorDynamics::Bus<ScalarT, IdxT> bus2(30.0, 40.0);
        nested.addBus(&bus2);

        PhasorDynamics::SignalNode<ScalarT, IdxT> bus1_vr;
        PhasorDynamics::SignalNode<ScalarT, IdxT> bus1_vi;
        PhasorDynamics::SignalNode<ScalarT, IdxT> bus2_vr;
        PhasorDynamics::SignalNode<ScalarT, IdxT> bus2_vi;
        bus1.getPorts().out[PhasorDynamics::BusSignalOutputs::vr].connect(&bus1_vr);
        bus1.getPorts().out[PhasorDynamics::BusSignalOutputs::vi].connect(&bus1_vi);
        bus2.getPorts().out[PhasorDynamics::BusSignalOutputs::vr].connect(&bus2_vr);
        bus2.getPorts().out[PhasorDynamics::BusSignalOutputs::vi].connect(&bus2_vi);

        // The branch lives in the nested system and connects buses across
        // system levels.
        PhasorDynamics::Branch<ScalarT, IdxT> branch(R, X, G, B);
        branch.getPorts().terminal_vr[PhasorDynamics::BranchBuses::bus1].connect(&bus1_vr);
        branch.getPorts().terminal_vi[PhasorDynamics::BranchBuses::bus1].connect(&bus1_vi);
        branch.getPorts().terminal_vr[PhasorDynamics::BranchBuses::bus2].connect(&bus2_vr);
        branch.getPorts().terminal_vi[PhasorDynamics::BranchBuses::bus2].connect(&bus2_vi);
        nested.addComponent(&branch);

        root.addComponent(&nested);

        root.allocate();
        root.initialize();
        root.evaluateResidual();

        success *= isEqual(branch.getExternalResidual()[0], Ir1);
        success *= isEqual(branch.getExternalResidual()[1], Ii1);
        success *= isEqual(bus2.Ir(), Ir2);
        success *= isEqual(bus2.Ii(), Ii2);

        // The nested system's variables live in root storage.
        success *= isEqual(root.y().getData()[0], 30.0);
        success *= isEqual(root.y().getData()[1], 40.0);

        // The nested system completes the rows it owns on its own.
        success *= (nested.evaluateResidual() == 0);
        success *= isEqual(bus2.Ir(), Ir2);
        success *= isEqual(bus2.Ii(), Ii2);

        return success.report(__func__);
      }

      /// A nested system at nonzero offset matches its flat twin
      TestOutcome nestedTwin()
      {
        TestStatus success = true;

        RealT R{2.0}; ///< Branch series resistance
        RealT X{4.0}; ///< Branch series reactance
        RealT G{0.2}; ///< Branch shunt conductance
        RealT B{1.2}; ///< Branch shunt charging

        // Both systems hold a slack bus and two PQ buses in the same
        // variable order; the twin nests the second PQ bus and its branch.
        auto build = [&](PhasorDynamics::SystemModel<ScalarT, IdxT>& head,
                         PhasorDynamics::SystemModel<ScalarT, IdxT>& tail,
                         PhasorDynamics::BusInfinite<ScalarT, IdxT>& bus1,
                         PhasorDynamics::Bus<ScalarT, IdxT>&         bus2,
                         PhasorDynamics::Bus<ScalarT, IdxT>&         bus3,
                         PhasorDynamics::Branch<ScalarT, IdxT>&      branch12,
                         PhasorDynamics::Branch<ScalarT, IdxT>&      branch23,
                         PhasorDynamics::SignalNode<ScalarT, IdxT>*  nodes)
        {
          head.addBus(&bus1);
          head.addBus(&bus2);

          bus1.getPorts().out[PhasorDynamics::BusSignalOutputs::vr].connect(&nodes[0]);
          bus1.getPorts().out[PhasorDynamics::BusSignalOutputs::vi].connect(&nodes[1]);
          bus2.getPorts().out[PhasorDynamics::BusSignalOutputs::vr].connect(&nodes[2]);
          bus2.getPorts().out[PhasorDynamics::BusSignalOutputs::vi].connect(&nodes[3]);
          bus3.getPorts().out[PhasorDynamics::BusSignalOutputs::vr].connect(&nodes[4]);
          bus3.getPorts().out[PhasorDynamics::BusSignalOutputs::vi].connect(&nodes[5]);

          branch12.getPorts().terminal_vr[PhasorDynamics::BranchBuses::bus1].connect(&nodes[0]);
          branch12.getPorts().terminal_vi[PhasorDynamics::BranchBuses::bus1].connect(&nodes[1]);
          branch12.getPorts().terminal_vr[PhasorDynamics::BranchBuses::bus2].connect(&nodes[2]);
          branch12.getPorts().terminal_vi[PhasorDynamics::BranchBuses::bus2].connect(&nodes[3]);
          head.addComponent(&branch12);

          branch23.getPorts().terminal_vr[PhasorDynamics::BranchBuses::bus1].connect(&nodes[2]);
          branch23.getPorts().terminal_vi[PhasorDynamics::BranchBuses::bus1].connect(&nodes[3]);
          branch23.getPorts().terminal_vr[PhasorDynamics::BranchBuses::bus2].connect(&nodes[4]);
          branch23.getPorts().terminal_vi[PhasorDynamics::BranchBuses::bus2].connect(&nodes[5]);
          tail.addBus(&bus3);
          tail.addComponent(&branch23);
        };

        // Flat reference: tail IS the head, so every component is direct.
        PhasorDynamics::SystemModel<ScalarT, IdxT> flat;
        PhasorDynamics::BusInfinite<ScalarT, IdxT> fbus1(10.0, 20.0);
        PhasorDynamics::Bus<ScalarT, IdxT>         fbus2(30.0, 40.0);
        PhasorDynamics::Bus<ScalarT, IdxT>         fbus3(25.0, 35.0);
        PhasorDynamics::Branch<ScalarT, IdxT>      fbranch12(R, X, G, B);
        PhasorDynamics::Branch<ScalarT, IdxT>      fbranch23(R, X, G, B);
        PhasorDynamics::SignalNode<ScalarT, IdxT>  fnodes[6];
        build(flat, flat, fbus1, fbus2, fbus3, fbranch12, fbranch23, fnodes);

        // Hierarchical twin: bus-3 and branch 2-3 nest at offset 2.
        PhasorDynamics::SystemModel<ScalarT, IdxT> root;
        PhasorDynamics::SystemModel<ScalarT, IdxT> nested;
        PhasorDynamics::BusInfinite<ScalarT, IdxT> nbus1(10.0, 20.0);
        PhasorDynamics::Bus<ScalarT, IdxT>         nbus2(30.0, 40.0);
        PhasorDynamics::Bus<ScalarT, IdxT>         nbus3(25.0, 35.0);
        PhasorDynamics::Branch<ScalarT, IdxT>      nbranch12(R, X, G, B);
        PhasorDynamics::Branch<ScalarT, IdxT>      nbranch23(R, X, G, B);
        PhasorDynamics::SignalNode<ScalarT, IdxT>  nnodes[6];
        build(root, nested, nbus1, nbus2, nbus3, nbranch12, nbranch23, nnodes);
        root.addComponent(&nested);

        flat.allocate();
        flat.initialize();
        flat.evaluateResidual();
        flat.tagDifferentiable();

        root.allocate();
        root.initialize();
        root.evaluateResidual();
        root.tagDifferentiable();

        success *= (root.size() == flat.size());

        // States, residuals, and tags match entry for entry.
        for (IdxT j = 0; j < flat.size(); ++j)
        {
          success *= isEqual(root.y().getData()[j], flat.y().getData()[j]);
          success *= isEqual(root.getResidual().getData()[j], flat.getResidual().getData()[j]);
          success *= (root.tag()[j] == flat.tag()[j]);
        }

        // The nested system's own tag vector is indexed locally.
        for (IdxT j = 0; j < nested.size(); ++j)
        {
          success *= (nested.tag()[j] == flat.tag()[2 + j]);
        }

#ifdef GRIDKIT_ENABLE_ENZYME
        // Jacobians assemble to the same CSR matrix. Allocation primed the
        // pattern, so these evaluations exercise the refill path.
        flat.evaluateJacobian();
        root.evaluateJacobian();

        auto* flat_jac = flat.getCsrJacobian();
        auto* root_jac = root.getCsrJacobian();

        success *= (flat_jac != nullptr && root_jac != nullptr);
        if (flat_jac != nullptr && root_jac != nullptr)
        {
          success *= (root_jac->getNnz() == flat_jac->getNnz());
          for (IdxT i = 0; i < flat.size() + 1; ++i)
          {
            success *= (root_jac->getRowData()[i] == flat_jac->getRowData()[i]);
          }
          for (IdxT i = 0; i < flat_jac->getNnz(); ++i)
          {
            success *= (root_jac->getColData()[i] == flat_jac->getColData()[i]);
            success *= isEqual(root_jac->getValues()[i], flat_jac->getValues()[i]);
          }
        }
#endif

        // Branch 2-3 injects into bus-2, which lies outside the nested
        // system, so that injection is the nested system's external residual.
        const auto& rows   = nested.getExternalResidualIndices();
        const auto  f_ext  = nested.getExternalResidual();
        const auto  f_23   = nbranch23.getExternalResidual();
        success           *= (rows.size() == 4);
        success           *= (rows[0] == nbus2.getResidualIndex(0));
        success           *= (rows[1] == nbus2.getResidualIndex(1));
        success           *= (rows[2] == INVALID_INDEX<IdxT>);
        success           *= (rows[3] == INVALID_INDEX<IdxT>);
        success           *= isEqual(f_ext[0], f_23[0]);
        success           *= isEqual(f_ext[1], f_23[1]);

        // A nested system evaluates its Jacobian entries; the root assembles.
        success *= (nested.evaluateJacobian() == 0);
        success *= (nested.getCsrJacobian() == nullptr);

        return success.report(__func__);
      }

      /// Systems nested two deep, next to an empty one, match their flat twin
      TestOutcome nestedDepthTwo()
      {
        TestStatus success = true;

        // A slack bus and three PQ buses in a chain, with the variables in
        // the same order in both systems. The twin nests bus-3 one level down
        // and bus-4 two levels down, each with the branch to the level above.
        struct Network
        {
          PhasorDynamics::BusInfinite<ScalarT, IdxT> bus1{10.0, 20.0};
          PhasorDynamics::Bus<ScalarT, IdxT>         bus2{30.0, 40.0};
          PhasorDynamics::Bus<ScalarT, IdxT>         bus3{25.0, 35.0};
          PhasorDynamics::Bus<ScalarT, IdxT>         bus4{15.0, 45.0};
          PhasorDynamics::Branch<ScalarT, IdxT>      branch12{2.0, 4.0, 0.2, 1.2};
          PhasorDynamics::Branch<ScalarT, IdxT>      branch23{2.0, 4.0, 0.2, 1.2};
          PhasorDynamics::Branch<ScalarT, IdxT>      branch34{2.0, 4.0, 0.2, 1.2};
          PhasorDynamics::SignalNode<ScalarT, IdxT>  nodes[8];

          void connect(PhasorDynamics::BusBase<ScalarT, IdxT>& bus, IdxT node)
          {
            bus.getPorts().out[PhasorDynamics::BusSignalOutputs::vr].connect(&nodes[node]);
            bus.getPorts().out[PhasorDynamics::BusSignalOutputs::vi].connect(&nodes[node + 1]);
          }

          void connect(PhasorDynamics::Branch<ScalarT, IdxT>& branch, BranchBuses terminal, IdxT node)
          {
            branch.getPorts().terminal_vr[terminal].connect(&nodes[node]);
            branch.getPorts().terminal_vi[terminal].connect(&nodes[node + 1]);
          }

          Network()
          {
            connect(bus1, 0);
            connect(bus2, 2);
            connect(bus3, 4);
            connect(bus4, 6);
            connect(branch12, BranchBuses::bus1, 0);
            connect(branch12, BranchBuses::bus2, 2);
            connect(branch23, BranchBuses::bus1, 2);
            connect(branch23, BranchBuses::bus2, 4);
            connect(branch34, BranchBuses::bus1, 4);
            connect(branch34, BranchBuses::bus2, 6);
          }
        };

        Network                                    f;
        PhasorDynamics::SystemModel<ScalarT, IdxT> flat;
        flat.addBus(&f.bus1);
        flat.addBus(&f.bus2);
        flat.addComponent(&f.branch12);
        flat.addBus(&f.bus3);
        flat.addComponent(&f.branch23);
        flat.addBus(&f.bus4);
        flat.addComponent(&f.branch34);

        Network                                    n;
        PhasorDynamics::SystemModel<ScalarT, IdxT> root;
        PhasorDynamics::SystemModel<ScalarT, IdxT> empty;
        PhasorDynamics::SystemModel<ScalarT, IdxT> middle;
        PhasorDynamics::SystemModel<ScalarT, IdxT> inner;
        inner.addBus(&n.bus4);
        inner.addComponent(&n.branch34);
        middle.addBus(&n.bus3);
        middle.addComponent(&n.branch23);
        middle.addComponent(&inner);
        root.addComponent(&empty);
        root.addBus(&n.bus1);
        root.addBus(&n.bus2);
        root.addComponent(&n.branch12);
        root.addComponent(&middle);

        flat.allocate();
        flat.initialize();
        flat.evaluateResidual();

        root.allocate();
        root.initialize();
        root.evaluateResidual();

        success *= (root.size() == flat.size());
        success *= (empty.size() == 0);
        for (IdxT j = 0; j < flat.size(); ++j)
        {
          success *= isEqual(root.y().getData()[j], flat.y().getData()[j]);
          success *= isEqual(root.getResidual().getData()[j], flat.getResidual().getData()[j]);
        }

        // Each subsystem passes up only its injections into the level above;
        // the root has nothing outside it.
        const auto& middle_rows  = middle.getExternalResidualIndices();
        const auto& inner_rows   = inner.getExternalResidualIndices();
        success                 *= (inner_rows[0] == n.bus3.getResidualIndex(0));
        success                 *= (inner_rows[2] == INVALID_INDEX<IdxT>);
        success                 *= (middle_rows[0] == n.bus2.getResidualIndex(0));
        success                 *= (middle_rows[2] == INVALID_INDEX<IdxT>);
        success                 *= (middle_rows[4] == INVALID_INDEX<IdxT>);
        for (const auto row : root.getExternalResidualIndices())
        {
          success *= (row == INVALID_INDEX<IdxT>);
        }

#ifdef GRIDKIT_ENABLE_ENZYME
        flat.evaluateJacobian();
        root.evaluateJacobian();

        auto* flat_jac = flat.getCsrJacobian();
        auto* root_jac = root.getCsrJacobian();

        success *= (flat_jac != nullptr && root_jac != nullptr);
        if (flat_jac != nullptr && root_jac != nullptr)
        {
          success *= (root_jac->getNnz() == flat_jac->getNnz());
          for (IdxT i = 0; i < flat_jac->getNnz(); ++i)
          {
            success *= (root_jac->getColData()[i] == flat_jac->getColData()[i]);
            success *= isEqual(root_jac->getValues()[i], flat_jac->getValues()[i]);
          }
        }
#endif

        return success.report(__func__);
      }

      /// A subsystem may read signals its parent produces after it
      TestOutcome nestedLaterProducer()
      {
        TestStatus success = true;

        using SignalIn  = PhasorDynamics::BusSignalInputs;
        using SignalOut = PhasorDynamics::BusSignalOutputs;

        // The nested bus takes its voltage from a bus the root adds after the
        // subsystem, and returns its boundary current to that bus.
        PhasorDynamics::SignalNode<ScalarT, IdxT> vr;
        PhasorDynamics::SignalNode<ScalarT, IdxT> vi;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ir;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ii;

        PhasorDynamics::BusSignalVoltageIn<ScalarT, IdxT>  bus_in;
        PhasorDynamics::BusSignalVoltageOut<ScalarT, IdxT> bus_out(1.02, -0.05);
        bus_in.getPorts().in.template port<SignalIn::vr>().connect(&vr);
        bus_in.getPorts().in.template port<SignalIn::vi>().connect(&vi);
        bus_in.getPorts().out.template port<SignalOut::ir>().connect(&ir);
        bus_in.getPorts().out.template port<SignalOut::ii>().connect(&ii);
        bus_out.getPorts().out.template port<SignalOut::vr>().connect(&vr);
        bus_out.getPorts().out.template port<SignalOut::vi>().connect(&vi);
        bus_out.getPorts().in.template port<SignalIn::ir>().connect(&ir);
        bus_out.getPorts().in.template port<SignalIn::ii>().connect(&ii);

        PhasorDynamics::SystemModel<ScalarT, IdxT> root;
        PhasorDynamics::SystemModel<ScalarT, IdxT> nested;
        nested.addBus(&bus_in);
        root.addComponent(&nested);
        root.addBus(&bus_out);

        // The root verifies once every producer in the tree is allocated.
        const auto previous_verbosity = Log::verbosity();
        // Suppress the error a failed verification would log.
        // Use EVERYTHING to inspect the diagnostic.
        Log::setVerbosity(Log::Verbosity::NONE);
        const bool allocated = !throws<std::runtime_error>(
            [&]()
            { root.allocate(); });
        Log::setVerbosity(previous_verbosity);

        success *= allocated;
        if (allocated)
        {
          success *= (root.evaluateResidual() == 0);
          success *= isEqual(bus_in.Vr(), 1.02);
          success *= isEqual(bus_in.Vi(), -0.05);
        }

        return success.report(__func__);
      }

      /// A subsystem adopts the system base of its parent
      TestOutcome nestedSystemBase()
      {
        TestStatus success = true;

        BaseProbeComponent                         probe;
        PhasorDynamics::SystemModel<ScalarT, IdxT> nested;
        nested.setSystemBase(50.0, 200.0e6);
        nested.addComponent(&probe);
        success *= isEqual(probe.vaSystemBase(), 200.0e6);

        PhasorDynamics::SystemModel<ScalarT, IdxT> root;
        root.setSystemBase(60.0, 100.0e6);
        root.addComponent(&nested);
        success *= isEqual(probe.vaSystemBase(), 100.0e6);

        return success.report(__func__);
      }

      TestOutcome reallocateAfterTopologyChange()
      {
        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT> system;
        PhasorDynamics::Bus<ScalarT, IdxT>         bus1(1.0, 0.0);
        PhasorDynamics::Bus<ScalarT, IdxT>         bus2(1.0, 0.0);
        PhasorDynamics::SignalNode<ScalarT, IdxT>  bus1_vr;
        PhasorDynamics::SignalNode<ScalarT, IdxT>  bus1_vi;
        PhasorDynamics::LoadZ<ScalarT, IdxT>       load;

        bus1.getPorts().out[PhasorDynamics::BusSignalOutputs::vr].connect(&bus1_vr);
        bus1.getPorts().out[PhasorDynamics::BusSignalOutputs::vi].connect(&bus1_vi);
        load.getPorts().terminal_vr[PhasorDynamics::LoadZBuses::bus].connect(&bus1_vr);
        load.getPorts().terminal_vi[PhasorDynamics::LoadZBuses::bus].connect(&bus1_vi);

        system.addBus(&bus1);
        system.addComponent(&load);
        success                    *= system.allocate() == 0;
        const IdxT size_before_bus  = system.size();

        system.addBus(&bus2);
        success *= system.allocate() == 0;
        success *= system.size() == size_before_bus + bus2.size();

#ifdef GRIDKIT_ENABLE_ENZYME
        const auto* jacobian  = system.getCsrJacobian();
        success              *= jacobian != nullptr;

        IdxT nnz_without_branch = 0;
        if (jacobian != nullptr)
        {
          success            *= jacobian->getNumRows() == system.size();
          success            *= jacobian->getNumColumns() == system.size();
          nnz_without_branch  = jacobian->getNnz();
        }
#endif

        PhasorDynamics::SignalNode<ScalarT, IdxT> bus2_vr;
        PhasorDynamics::SignalNode<ScalarT, IdxT> bus2_vi;
        bus2.getPorts().out[PhasorDynamics::BusSignalOutputs::vr].connect(&bus2_vr);
        bus2.getPorts().out[PhasorDynamics::BusSignalOutputs::vi].connect(&bus2_vi);

        PhasorDynamics::Branch<ScalarT, IdxT> branch;
        branch.getPorts().terminal_vr[PhasorDynamics::BranchBuses::bus1].connect(&bus1_vr);
        branch.getPorts().terminal_vi[PhasorDynamics::BranchBuses::bus1].connect(&bus1_vi);
        branch.getPorts().terminal_vr[PhasorDynamics::BranchBuses::bus2].connect(&bus2_vr);
        branch.getPorts().terminal_vi[PhasorDynamics::BranchBuses::bus2].connect(&bus2_vi);
        system.addComponent(&branch);
        success *= system.allocate() == 0;
        success *= system.evaluateJacobian() == 0;

#ifdef GRIDKIT_ENABLE_ENZYME
        jacobian  = system.getCsrJacobian();
        success  *= jacobian != nullptr;
        if (jacobian != nullptr)
        {
          success *= jacobian->getNumRows() == system.size();
          success *= jacobian->getNumColumns() == system.size();
          success *= jacobian->getNnz() > nnz_without_branch;
        }
#endif

        return success.report(__func__);
      }

      TestOutcome modelVectorsAliasSystemStorage()
      {
        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT> system;
        PhasorDynamics::Bus<ScalarT, IdxT>         bus1(1.0, 2.0);
        PhasorDynamics::BusInfinite<ScalarT, IdxT> infinite_bus;
        PhasorDynamics::Bus<ScalarT, IdxT>         bus2(3.0, 4.0);
        PhasorDynamics::SignalNode<ScalarT, IdxT>  bus2_vr;
        PhasorDynamics::SignalNode<ScalarT, IdxT>  bus2_vi;
        PhasorDynamics::SignalNode<ScalarT, IdxT>  bus1_vr;
        PhasorDynamics::SignalNode<ScalarT, IdxT>  bus1_vi;
        PhasorDynamics::Branch<ScalarT, IdxT>      branch;
        PhasorDynamics::LoadZ<ScalarT, IdxT>       load(1.0, 1.0);

        bus1.getPorts().out[PhasorDynamics::BusSignalOutputs::vr].connect(&bus1_vr);
        bus1.getPorts().out[PhasorDynamics::BusSignalOutputs::vi].connect(&bus1_vi);
        bus2.getPorts().out[PhasorDynamics::BusSignalOutputs::vr].connect(&bus2_vr);
        bus2.getPorts().out[PhasorDynamics::BusSignalOutputs::vi].connect(&bus2_vi);
        branch.getPorts().terminal_vr[PhasorDynamics::BranchBuses::bus1].connect(&bus1_vr);
        branch.getPorts().terminal_vi[PhasorDynamics::BranchBuses::bus1].connect(&bus1_vi);
        branch.getPorts().terminal_vr[PhasorDynamics::BranchBuses::bus2].connect(&bus2_vr);
        branch.getPorts().terminal_vi[PhasorDynamics::BranchBuses::bus2].connect(&bus2_vi);
        load.getPorts().terminal_vr[PhasorDynamics::LoadZBuses::bus].connect(&bus2_vr);
        load.getPorts().terminal_vi[PhasorDynamics::LoadZBuses::bus].connect(&bus2_vi);

        system.addBus(&bus1);
        system.addBus(&infinite_bus);
        system.addBus(&bus2);
        system.addComponent(&branch);
        system.addComponent(&load);

        if (system.allocate() != 0
            || system.setAbsoluteTolerance(1e-4) != 0)
        {
          success = false;
          return success.report(__func__);
        }

        auto checkAlias = [&](auto& system_vector, auto& model_vector, IdxT offset)
        {
          auto*      system_data = system_vector.getData();
          auto*      model_data  = model_vector.getData();
          const auto first       = static_cast<std::size_t>(offset);

          if (!system_data || model_data != system_data + first)
          {
            success = false;
            return;
          }

          success *= system_vector.setToConst(ScalarT{3.0}) == 0;
          success *= isEqual(model_data[0], ScalarT{3.0});

          success *= model_vector.setToConst(ScalarT{4.0}) == 0;
          success *= isEqual(system_data[first], ScalarT{4.0});
        };

        auto checkModel = [&](auto& model, IdxT offset)
        {
          success *= model.getVariableIndex(0) == offset;
          success *= model.getResidualIndex(0) == offset;

          checkAlias(system.y(), model.y(), offset);
          checkAlias(system.yp(), model.yp(), offset);
          checkAlias(system.getResidual(), model.getResidual(), offset);
          checkAlias(system.absoluteTolerance(), model.absoluteTolerance(), offset);
        };

        const IdxT bus2_offset = bus1.size();
        const IdxT load_offset = bus1.size() + bus2.size();
        const auto bus2_first  = static_cast<std::size_t>(bus2_offset);

        auto rebind = [&](auto& model, IdxT offset)
        {
          return model.bind(system.y(),
                            system.yp(),
                            system.getResidual(),
                            system.absoluteTolerance(),
                            offset);
        };

        // Rebinding the same slices is a no-op.
        success *= rebind(bus2, bus2_offset) == 0;
        success *= rebind(load, load_offset) == 0;

        checkModel(bus2, bus2_offset);
        checkModel(load, load_offset);

        // Tags remain model-owned and are collected separately.
        system.tag()[bus2_first]  = true;
        success                  *= system.tagDifferentiable() == 0;
        success                  *= !system.tag()[bus2_first];

        bus2.tag()[0]  = true;
        success       *= !system.tag()[bus2_first];

        return success.report(__func__);
      }

      /**
       * @brief Test for exception when signals are incorrectly configured
       */
      TestOutcome signalError()
      {
        using namespace GridKit::PhasorDynamics;
        auto data = parseSystemModelData("ThreeBusBasicBad.json");
        auto sys  = SystemModel<double, size_t>(data);

        TestStatus status{true};
        const auto previous_verbosity = Log::verbosity();
        // Suppress the expected signal-configuration error below.
        // Use EVERYTHING to inspect the diagnostic.
        Log::setVerbosity(Log::Verbosity::NONE);
        status *= throws<std::runtime_error>(
            [&]()
            { sys.allocate(); });
        Log::setVerbosity(previous_verbosity);

        return status.report(__func__);
      }

      /**
       * @brief Test for exception when a child cannot bind to system storage
       */
      TestOutcome allocationError()
      {
        using namespace GridKit::PhasorDynamics;

        TestStatus                 status{true};
        SystemModel<ScalarT, IdxT> system;
        Bus<ScalarT, IdxT>         bus(ScalarT{1.0}, ScalarT{0.0});

        status *= bus.allocate() == 0;
        system.addBus(&bus);
        const auto previous_verbosity = Log::verbosity();
        // Suppress the expected child-allocation error below.
        // Use EVERYTHING to inspect the diagnostic.
        Log::setVerbosity(Log::Verbosity::NONE);
        status *= throws<std::runtime_error>(
            [&]()
            { system.allocate(); });
        Log::setVerbosity(previous_verbosity);

        return status.report(__func__);
      }

      /// SystemModel propagates a statically valid component's initialization error.
      TestOutcome componentInitializationError()
      {
        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT> system;
        InitializationFailureComponent             component;
        system.addComponent(&component);

        success *= system.verify() == 0;

        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::NONE);

        if (system.hasJacobian())
        {
          success *= throws<std::runtime_error>([&]()
                                                { system.allocate(); });
        }
        else
        {
          success *= system.allocate() == 0;
          success *= system.initialize() != 0;
        }

        Log::setVerbosity(previous_verbosity);
        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      TestOutcome jacobian()
      {
        TestStatus success = true;

        PhasorDynamics::SystemModelData<ScalarT, IdxT> data;

        // Set bus data
        data.bus.resize(2);

        // Bus 0
        data.bus[0].bus_id   = 0;
        data.bus[0].bus_type = PhasorDynamics::BusData<ScalarT, IdxT>::BusType::SLACK;
        data.bus[0].Vr0      = 10.0;
        data.bus[0].Vi0      = 20.0;

        // Bus 1
        data.bus[1].bus_id   = 1;
        data.bus[1].bus_type = PhasorDynamics::BusData<ScalarT, IdxT>::BusType::DEFAULT;
        data.bus[1].Vr0      = 30.0;
        data.bus[1].Vi0      = 40.0;

        // Set branch data
        data.branch.resize(1);

        // Branch 0-1
        data.branch[0].buses[BranchBuses::bus1]        = data.bus[0].bus_id;
        data.branch[0].buses[BranchBuses::bus2]        = data.bus[1].bus_id;
        data.branch[0].parameters[BranchParameters::R] = 2.0;
        data.branch[0].parameters[BranchParameters::X] = 4.0;
        data.branch[0].parameters[BranchParameters::G] = 0.2;
        data.branch[0].parameters[BranchParameters::B] = 1.2;

        // Jacobian via DependencyTracking
        std::vector<DependencyTracking::Variable::DependencyMap> dependency_tracking_jacobian = DependencyTrackingJacobian(data);

        // Jacobian via Enzyme
        std::vector<DependencyTracking::Variable::DependencyMap> enzyme_jacobian = EnzymeJacobian(data);

        /// Compare DependencyTracking dependencies to Enzyme's
        for (size_t i = 0; i < dependency_tracking_jacobian.size(); ++i)
        {
          success *= (GridKit::Testing::isEqual(dependency_tracking_jacobian[i], enzyme_jacobian[i]));
        }

        return success.report(__func__);
      }

    private:
      std::vector<DependencyTracking::Variable::DependencyMap> DependencyTrackingJacobian(
          PhasorDynamics::SystemModelData<ScalarT, IdxT> data)
      {
        // Create an empty system model
        PhasorDynamics::SystemModel<DependencyTracking::Variable, IdxT> system(data);

        // Allocate and initialize the system
        system.allocate();
        system.initialize();

        // Evaluate and get the system Jacobian
        system.evaluateResidual();
        system.evaluateJacobian();
        auto* system_jacobian = system.getCsrJacobian();
        auto& output_stream   = Log::misc();
        output_stream << "Sparse Csr Matrix: System Jacobian with DependencyTracking\n";
        system_jacobian->print(output_stream);

        return GridKit::Testing::MapFromCsr(system_jacobian);
      }

      std::vector<DependencyTracking::Variable::DependencyMap> EnzymeJacobian(
          PhasorDynamics::SystemModelData<ScalarT, IdxT> data)
      {
        // Create an empty system model
        PhasorDynamics::SystemModel<ScalarT, IdxT> system(data);

        // Allocate and initialize the system
        system.allocate();
        system.initialize();

        // Evaluate and get the system Jacobian
        system.evaluateResidual();
        system.evaluateJacobian();
        auto* system_jacobian = system.getCsrJacobian();
        auto& output_stream   = Log::misc();
        output_stream << "Sparse Csr Matrix: System Jacobian with Enzyme\n";
        system_jacobian->print(output_stream);

        return GridKit::Testing::MapFromCsr(system_jacobian);
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
