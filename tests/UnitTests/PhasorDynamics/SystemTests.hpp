#pragma once

#include <array>
#include <cstddef>
#include <cstring>
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
#include <GridKit/Model/PhasorDynamics/BusFault/BusFault.hpp>
#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/Load/LoadZ/LoadZ.hpp>
#include <GridKit/Model/PhasorDynamics/Load/LoadZIP/LoadZIP.hpp>
#include <GridKit/Model/PhasorDynamics/NetworkAdmittance.hpp>
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

      class EvaluationContextProbe final : public PhasorDynamics::LoadZ<ScalarT, IdxT>
      {
        using BaseT = PhasorDynamics::LoadZ<ScalarT, IdxT>;

      public:
        using BaseT::BaseT;

        const RealT& evaluationTime() const
        {
          return this->time();
        }

        const RealT& evaluationAlpha() const
        {
          return this->alpha();
        }

        std::uint64_t evaluationAdmittanceEpoch() const
        {
          return this->admittanceEpoch();
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

        auto* bus0 = system->getBus(0);
        auto* bus1 = system->getBus(1);

        success *= isEqual(bus0->Ir(), Ir0);
        success *= isEqual(bus0->Ii(), Ii0);
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

        PhasorDynamics::Branch<ScalarT, IdxT> branch(&bus1, &bus2, R, X, G, B);
        system.addComponent(&branch);

        system.allocate();
        system.initialize();
        system.evaluateResidual();

        success *= isEqual(bus1.Ir(), Ir1);
        success *= isEqual(bus1.Ii(), Ii1);
        success *= isEqual(bus2.Ir(), Ir2);
        success *= isEqual(bus2.Ii(), Ii2);

        return success.report(__func__);
      }

      TestOutcome residualAssemblyIsIdempotent()
      {
        using namespace PhasorDynamics;

        TestStatus success = true;

        SystemModel<ScalarT, IdxT> system;
        Bus<ScalarT, IdxT>         bus1(ScalarT{10.0}, ScalarT{20.0});
        Bus<ScalarT, IdxT>         bus2(ScalarT{30.0}, ScalarT{40.0});
        Branch<ScalarT, IdxT>      branch(&bus1,
                                     &bus2,
                                     RealT{2.0},
                                     RealT{4.0},
                                     RealT{0.2},
                                     RealT{1.2});

        bus1.setBusID(IdxT{0});
        bus2.setBusID(IdxT{1});
        system.addBus(&bus1);
        system.addBus(&bus2);
        system.addComponent(&branch);

        if (system.allocate() != 0 || system.initialize() != 0)
        {
          success = false;
          return success.report(__func__);
        }

        const std::array<ScalarT, 4> expected{
            ScalarT{17.0},
            ScalarT{-10.0},
            ScalarT{15.0},
            ScalarT{-20.0}};

        if (system.evaluateResidual() != 0)
        {
          success = false;
          return success.report(__func__);
        }

        const auto* first_data = system.getResidual().getData(memory::HOST);
        if (first_data == nullptr)
        {
          success = false;
          return success.report(__func__);
        }

        std::array<ScalarT, 4> first{};
        std::memcpy(first.data(), first_data, first.size() * sizeof(ScalarT));

        if (system.evaluateResidual() != 0)
        {
          success = false;
          return success.report(__func__);
        }

        const auto* second_data = system.getResidual().getData(memory::HOST);
        if (second_data == nullptr)
        {
          success = false;
          return success.report(__func__);
        }

        success *= std::memcmp(first.data(),
                               second_data,
                               first.size() * sizeof(ScalarT))
                   == 0;
        for (std::size_t i = 0; i < expected.size(); ++i)
        {
          success *= isEqual(second_data[i], expected[i]);
        }

        // System residual assembly restores HOST freshness at both alias levels.
        system.getResidual().setDataUpdated(memory::DEVICE);
        bus1.getResidual().setDataUpdated(memory::DEVICE);
        bus2.getResidual().setDataUpdated(memory::DEVICE);

        success *= system.evaluateResidual() == 0;
        success *= system.getResidual().getData(memory::HOST) != nullptr;
        success *= bus1.getResidual().getData(memory::HOST) != nullptr;
        success *= bus2.getResidual().getData(memory::HOST) != nullptr;

        return success.report(__func__);
      }

      TestOutcome networkAdmittanceMergesDuplicateStamps()
      {
        using namespace PhasorDynamics;

        TestStatus success = true;
        using StampT       = AdmittanceStamp<RealT, IdxT>;

        NetworkAdmittance<ScalarT, IdxT> network;
        std::vector<StampT> stamps{
            {IdxT{0}, IdxT{0}, RealT{1.0}, RealT{2.0}},
            {IdxT{0}, IdxT{0}, RealT{3.0}, RealT{4.0}},
            {IdxT{0}, IdxT{2}, RealT{-1.0}, RealT{0.5}},
            {IdxT{2}, IdxT{0}, RealT{2.0}, RealT{-1.0}}};

        network.assemble(stamps, {IdxT{0}, IdxT{2}, IdxT{4}});
        success *= network.rowCount() == IdxT{3};
        success *= network.nnz() == IdxT{3};

        const std::array<ScalarT, 6> y{
            ScalarT{1.0},
            ScalarT{2.0},
            ScalarT{3.0},
            ScalarT{4.0},
            ScalarT{5.0},
            ScalarT{6.0}};
        std::array<ScalarT, 6> f{
            ScalarT{9.0},
            ScalarT{9.0},
            ScalarT{9.0},
            ScalarT{9.0},
            ScalarT{9.0},
            ScalarT{9.0}};

        network.multiply(y.data(), f.data());
        const std::array<ScalarT, 6> expected{
            ScalarT{-13.0},
            ScalarT{11.5},
            ScalarT{4.0},
            ScalarT{3.0},
            ScalarT{0.0},
            ScalarT{0.0}};
        for (std::size_t i = 0; i < expected.size(); ++i)
        {
          success *= isEqual(f[i], expected[i]);
        }

        NetworkAdmittance<ScalarT, IdxT> invalid_network;
        std::vector<StampT> bad_stamps{
            {IdxT{1}, IdxT{0}, RealT{1.0}, RealT{0.0}}};
        success *= throws<std::invalid_argument>(
            [&]()
            { invalid_network.assemble(bad_stamps, {IdxT{0}}); });

        return success.report(__func__);
      }

      TestOutcome networkAdmittanceTracksParameterChanges()
      {
        using namespace PhasorDynamics;

        TestStatus success = true;

        SystemModel<ScalarT, IdxT> system;
        Bus<ScalarT, IdxT>         bus(ScalarT{1.0}, ScalarT{0.0});
        LoadZ<ScalarT, IdxT>       load_z(&bus, RealT{1.0}, RealT{0.0});
        LoadZIP<ScalarT, IdxT>     load_zip(
            &bus, RealT{2.0}, RealT{1.0}, RealT{0.0}, RealT{0.0});

        bus.setBusID(IdxT{0});
        system.addBus(&bus);
        system.addComponent(&load_z);
        system.addComponent(&load_zip);
        if (system.allocate() != 0 || system.initialize() != 0)
        {
          success = false;
          return success.report(__func__);
        }

        auto expectBusResidual = [&](ScalarT expected_ir, ScalarT expected_ii)
        {
          success *= system.evaluateResidual() == 0;
          success *= isEqual(bus.Ir(), expected_ir);
          success *= isEqual(bus.Ii(), expected_ii);
        };

        // LoadZ and pure-Z LoadZIP stamp the same entry and are merged.
        expectBusResidual(ScalarT{-3.0}, ScalarT{1.0});

        // Rebuilds overwrite stale HOST residual output instead of reading it.
        system.getResidual().setDataUpdated(memory::DEVICE);
        load_z.setR(RealT{2.0});
        load_zip.setPnom(RealT{4.0});
        expectBusResidual(ScalarT{-4.5}, ScalarT{1.0});

        // A mixed ZIP leaves the matrix and returns to the component sweep.
        load_zip.setAlphaI(RealT{0.25});
        expectBusResidual(ScalarT{-4.5}, ScalarT{1.0});

        bus.Vr() = ScalarT{2.0};
        bus.Vi() = ScalarT{0.0};
        system.y().setDataUpdated(memory::HOST);
        expectBusResidual(ScalarT{-8.0}, ScalarT{1.75});

        // Returning to pure-Z stamps the load again with current parameters.
        load_zip.setAlphaI(RealT{0.0});
        expectBusResidual(ScalarT{-9.0}, ScalarT{2.0});

        // Reinitialization also reassembles without reading stale residuals.
        system.getResidual().setDataUpdated(memory::DEVICE);
        success *= system.initialize() == 0;
        expectBusResidual(ScalarT{-4.5}, ScalarT{1.0});

        SystemModel<ScalarT, IdxT> branch_system;
        Bus<ScalarT, IdxT>         bus1(ScalarT{1.0}, ScalarT{0.0});
        Bus<ScalarT, IdxT>         bus2(ScalarT{0.0}, ScalarT{0.0});
        Branch<ScalarT, IdxT>      branch(
            &bus1, &bus2, RealT{1.0}, RealT{0.0}, RealT{0.0}, RealT{0.0});

        bus1.setBusID(IdxT{1});
        bus2.setBusID(IdxT{2});
        branch_system.addBus(&bus1);
        branch_system.addBus(&bus2);
        branch_system.addComponent(&branch);
        success *= branch_system.allocate() == 0;
        success *= branch_system.initialize() == 0;
        success *= branch_system.evaluateResidual() == 0;
        success *= isEqual(bus1.Ir(), ScalarT{-1.0});
        success *= isEqual(bus2.Ir(), ScalarT{1.0});

        branch.setR(RealT{2.0});
        success *= branch_system.evaluateResidual() == 0;
        success *= isEqual(bus1.Ir(), ScalarT{-0.5});
        success *= isEqual(bus2.Ir(), ScalarT{0.5});

        return success.report(__func__);
      }

      TestOutcome networkAdmittancePreservesFaultEvents()
      {
        using namespace PhasorDynamics;

        TestStatus success = true;

        SystemModel<ScalarT, IdxT> system;
        Bus<ScalarT, IdxT>         bus(ScalarT{1.0}, ScalarT{0.0});
        LoadZ<ScalarT, IdxT>       load(&bus, RealT{1.0}, RealT{0.0});
        BusFault<ScalarT, IdxT>    fault(&bus, RealT{1.0}, RealT{0.0}, 0);

        bus.setBusID(IdxT{0});
        system.addBus(&bus);
        system.addComponent(&load);
        system.addFault(&fault);
        if (system.allocate() != 0 || system.initialize() != 0)
        {
          success = false;
          return success.report(__func__);
        }

        success *= system.evaluateResidual() == 0;
        success *= isEqual(bus.Ir(), ScalarT{-1.0});
        success *= isEqual(bus.Ii(), ScalarT{0.0});

        auto* fault_state = fault.y().getData(memory::HOST);
        if (fault_state == nullptr)
        {
          success = false;
          return success.report(__func__);
        }
        fault_state[0] = ScalarT{-1.0};
        fault_state[1] = ScalarT{0.0};
        fault.y().setDataUpdated(memory::HOST);

        fault.setStatus(true);
        success *= system.evaluateResidual() == 0;
        success *= isEqual(bus.Ir(), ScalarT{-2.0});
        success *= isEqual(bus.Ii(), ScalarT{0.0});

        fault.setStatus(false);
        success *= system.evaluateResidual() == 0;
        success *= isEqual(bus.Ir(), ScalarT{-1.0});
        success *= isEqual(bus.Ii(), ScalarT{0.0});

        return success.report(__func__);
      }

      TestOutcome reallocateAfterTopologyChange()
      {
        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT> system;
        PhasorDynamics::Bus<ScalarT, IdxT>         bus1(1.0, 0.0);
        PhasorDynamics::Bus<ScalarT, IdxT>         bus2(1.0, 0.0);
        PhasorDynamics::BusFault<ScalarT, IdxT>    fault(&bus1);

        system.addBus(&bus1);
        system.addComponent(&fault);
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

        PhasorDynamics::Branch<ScalarT, IdxT> branch(&bus1, &bus2);
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
        PhasorDynamics::Branch<ScalarT, IdxT>      branch(&bus1, &bus2);
        PhasorDynamics::BusFault<ScalarT, IdxT>    fault(&bus2);

        system.addBus(&bus1);
        system.addBus(&infinite_bus);
        system.addBus(&bus2);
        system.addComponent(&branch);
        system.addComponent(&fault);

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

        const IdxT bus2_offset  = bus1.size();
        const IdxT fault_offset = bus1.size() + bus2.size();
        const auto bus2_first   = static_cast<std::size_t>(bus2_offset);

        auto rebind = [&](auto& model, IdxT offset)
        {
          return model.bind(system.y(),
                            system.yp(),
                            system.getResidual(),
                            system.absoluteTolerance(),
                            offset);
        };

        // Rebinding to a different slice refreshes the bus terminal aliases.
        auto* const original_Vr  = &bus2.Vr();
        success                 *= rebind(bus2, IdxT{0}) == 0;
        success                 *= &bus2.Vr() == system.y().getData();
        success                 *= &bus2.Ir() == system.getResidual().getData();
        success                 *= &bus2.Vr() != original_Vr;

        success *= rebind(bus2, bus2_offset) == 0;
        success *= &bus2.Vr() == system.y().getData() + bus2_offset;
        success *= &bus2.Ir() == system.getResidual().getData() + bus2_offset;

        // Rebinding the remaining model to the same slice is a no-op.
        success *= rebind(fault, fault_offset) == 0;

        checkModel(bus2, bus2_offset);
        checkModel(fault, fault_offset);

        // Tags remain model-owned and are collected separately.
        system.tag()[bus2_first]  = true;
        success                  *= system.tagDifferentiable() == 0;
        success                  *= !system.tag()[bus2_first];

        bus2.tag()[0]  = true;
        success       *= !system.tag()[bus2_first];

        return success.report(__func__);
      }

      TestOutcome componentsShareEvaluationContext()
      {
        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT> system;
        PhasorDynamics::Bus<ScalarT, IdxT>         bus(1.0, 0.0);
        EvaluationContextProbe                     load1(&bus, 1.0, 1.0);
        EvaluationContextProbe                     load2(&bus, 1.0, 1.0);

        load1.updateTime(0.25, 0.5);
        success *= isEqual(load1.evaluationTime(), RealT{0.25});
        success *= isEqual(load1.evaluationAlpha(), RealT{0.5});

        system.addBus(&bus);
        system.addComponent(&load1);
        system.addComponent(&load2);
        success *= system.allocate() == 0;

        system.updateTime(1.25, 2.5);
        success *= isEqual(load1.evaluationTime(), RealT{1.25});
        success *= isEqual(load1.evaluationAlpha(), RealT{2.5});
        success *= isEqual(load2.evaluationTime(), RealT{1.25});
        success *= isEqual(load2.evaluationAlpha(), RealT{2.5});

        const auto epoch = load2.evaluationAdmittanceEpoch();
        load1.setR(2.0);
        success *= load1.evaluationAdmittanceEpoch() == epoch + 1;
        success *= load2.evaluationAdmittanceEpoch() == epoch + 1;

        system.updateTime(3.0, 4.0);
        success *= isEqual(load1.evaluationTime(), RealT{3.0});
        success *= isEqual(load1.evaluationAlpha(), RealT{4.0});
        success *= isEqual(load2.evaluationTime(), RealT{3.0});
        success *= isEqual(load2.evaluationAlpha(), RealT{4.0});

        return success.report(__func__);
      }


      /**
       * @brief Test for exception when signals are incorrectly configured
       */
      TestOutcome signalError()
      {
        using namespace std::filesystem;
        using namespace GridKit::PhasorDynamics;
        auto input_file = current_path() / "ThreeBusBasicBad.json";
        auto data       = parseSystemModelData(input_file);
        auto sys        = SystemModel<double, size_t>(data);

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
      TestOutcome jacobianAssemblyTracksCachedContributions()
      {
        using namespace PhasorDynamics;

        TestStatus success = true;

        SystemModel<ScalarT, IdxT> system;
        Bus<ScalarT, IdxT>         bus(ScalarT{1.0}, ScalarT{0.0});
        LoadZ<ScalarT, IdxT>       load(&bus, RealT{1.0}, RealT{0.0});
        LoadZIP<ScalarT, IdxT>     zip_load(
            &bus, RealT{2.0}, RealT{1.0}, RealT{0.0}, RealT{0.0});
        BusFault<ScalarT, IdxT>    fault(&bus, RealT{1.0}, RealT{0.0}, 0);

        bus.setBusID(IdxT{0});
        system.addBus(&bus);
        system.addComponent(&load);
        system.addComponent(&zip_load);
        system.addFault(&fault);
        if (system.allocate() != 0 || system.initialize() != 0)
        {
          success = false;
          return success.report(__func__);
        }

        auto valueAt = [&](IdxT row, IdxT col, RealT& value)
        {
          auto* jacobian = system.getCsrJacobian();
          if (jacobian == nullptr || row >= jacobian->getNumRows())
          {
            return false;
          }

          const auto* row_ptrs = jacobian->getRowData();
          const auto* columns  = jacobian->getColData();
          const auto* values   = jacobian->getValues();
          for (IdxT i = row_ptrs[row]; i < row_ptrs[row + 1]; ++i)
          {
            if (columns[i] == col)
            {
              value = values[i];
              return true;
            }
          }
          return false;
        };

        auto expectValue = [&](IdxT row, IdxT col, RealT expected)
        {
          RealT value{0};
          const bool found = valueAt(row, col, value);
          success         *= found;
          if (found)
          {
            success *= isEqual(value, expected);
          }
        };

        const IdxT bus_ir         = bus.getResidualIndex(0);
        const IdxT bus_vr         = bus.getVariableIndex(0);
        const IdxT fault_current_r = fault.getVariableIndex(0);

        // The initial load is invariant while the inactive fault contribution
        // remains a structural zero in the varying map.
        success *= system.evaluateJacobian() == 0;
        expectValue(bus_ir, bus_vr, RealT{-3.0});
        expectValue(bus_ir, fault_current_r, RealT{0.0});

        auto* first_jacobian = system.getCsrJacobian();
        if (first_jacobian == nullptr)
        {
          success = false;
          return success.report(__func__);
        }
        std::vector<RealT> first_values(
            first_jacobian->getValues(),
            first_jacobian->getValues() + first_jacobian->getNnz());

        success *= system.evaluateJacobian() == 0;
        auto* repeated_jacobian = system.getCsrJacobian();
        success *= repeated_jacobian->getNnz()
                   == static_cast<IdxT>(first_values.size());
        for (IdxT i = 0; i < repeated_jacobian->getNnz(); ++i)
        {
          success *= isEqual(repeated_jacobian->getValues()[i], first_values[i]);
        }

        // Dynamic blocks are still evaluated on every call.
        fault.setStatus(true);
        success *= system.evaluateJacobian() == 0;
        expectValue(bus_ir, bus_vr, RealT{-3.0});
        expectValue(bus_ir, fault_current_r, RealT{1.0});

        // An admittance mutation rebuilds the network and resnapshots the
        // constant Jacobian without losing the current dynamic contribution.
        load.setR(RealT{2.0});
        success *= system.evaluateJacobian() == 0;
        expectValue(bus_ir, bus_vr, RealT{-2.5});
        expectValue(bus_ir, fault_current_r, RealT{1.0});

        fault.setStatus(false);
        success *= system.evaluateJacobian() == 0;
        expectValue(bus_ir, bus_vr, RealT{-2.5});
        expectValue(bus_ir, fault_current_r, RealT{0.0});

        // A pure-Z ZIP moving into the dynamic sweep is immediately reflected
        // in the baseline/source map, then follows voltage changes per call.
        zip_load.setAlphaP(RealT{0.25});
        success *= system.evaluateJacobian() == 0;
        expectValue(bus_ir, bus_vr, RealT{-1.5});

        bus.Vr() = ScalarT{2.0};
        system.y().setDataUpdated(memory::HOST);
        success *= system.evaluateJacobian() == 0;
        expectValue(bus_ir, bus_vr, RealT{-1.875});

        // Returning to pure Z moves the ZIP contribution back into the
        // invariant snapshot.
        zip_load.setAlphaP(RealT{0.0});
        success *= system.evaluateJacobian() == 0;
        expectValue(bus_ir, bus_vr, RealT{-2.5});

        // Reallocation clears all structure, pointer, and baseline caches even
        // when a new zero-state stamp does not change the system dimension.
        LoadZ<ScalarT, IdxT> extra_load(&bus, RealT{2.0}, RealT{0.0});
        system.addComponent(&extra_load);
        success *= system.allocate() == 0;
        success *= system.initialize() == 0;
        success *= system.evaluateJacobian() == 0;
        expectValue(bus_ir, bus_vr, RealT{-3.0});
        expectValue(bus_ir, fault_current_r, RealT{0.0});

        return success.report(__func__);
      }

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
        std::cout << "Sparse Csr Matrix: System Jacobian with DependencyTracking\n";
        system_jacobian->print();

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
        std::cout << "Sparse Csr Matrix: System Jacobian with Enzyme\n";
        system_jacobian->print();

        return GridKit::Testing::MapFromCsr(system_jacobian);
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
