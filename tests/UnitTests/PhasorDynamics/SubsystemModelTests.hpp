#pragma once

#include <exception>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelDataJSONParser.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    using json = nlohmann::json;
    using Log  = ::GridKit::Utilities::Logger;

    template <class ScalarT, typename IdxT>
    class SubsystemModelTests
    {
    private:
      using ComponentT = PhasorDynamics::Component<ScalarT, IdxT>;
      using RealT      = typename ComponentT::RealT;
      using DataT      = PhasorDynamics::SystemModelData<RealT, IdxT>;
      using SystemT    = PhasorDynamics::SystemModel<ScalarT, IdxT>;
      using BusT       = PhasorDynamics::BusBase<ScalarT, IdxT>;

    public:
      SubsystemModelTests()  = default;
      ~SubsystemModelTests() = default;

      /**
       * @brief Compare a flat system with the same system split into
       * hierarchical subsystems built programmatically.
       */
      TestOutcome flatVsHierarchical()
      {
        TestStatus success = true;

        // Data of each system of the hierarchical input without its
        // subsystems: top, A, B, C
        const json         input   = readJson("SubsystemHierarchical.json");
        const json&        b_input = input["subsystems"][1];
        std::vector<DataT> data;
        for (json system : {input, input["subsystems"][0], b_input, b_input["subsystems"][0]})
        {
          system.erase("subsystems");
          system["header"] = input["header"];
          system["params"] = input["params"];
          data.push_back(system.template get<DataT>());
        }

        SystemT flat(readJson("SubsystemFlat.json").template get<DataT>());

        // The top system and B own (and delete) their subsystems. B contains
        // only bus 3 and C, so C gets branch 2 to bus 1.
        SystemT top(data[0]);
        auto*   a = new SystemT(data[1], {top.getBus(1)});
        auto*   b = new SystemT(data[2], {top.getBus(1)});
        auto*   c = new SystemT(data[3], {top.getBus(1), b->getBus(3)});
        b->addSubsystem(c);
        top.addSubsystem(a);
        top.addSubsystem(b);

        compareSystems(flat, top, success);

        return success.report(__func__);
      }

      /**
       * @brief Build the hierarchical system of flatVsHierarchical() from
       * hierarchical input data.
       */
      TestOutcome hierarchicalInput()
      {
        TestStatus success = true;

        const DataT top_data = readJson("SubsystemHierarchical.json").template get<DataT>();

        success *= top_data.subsystems.size() == 2;
        success *= top_data.subsystems[1].subsystems.size() == 1;
        success *= top_data.subsystems[1].subsystems[0].va_base == top_data.va_base;
        success *= top_data.level == 0;
        success *= top_data.subsystems[1].id == "B";
        success *= top_data.subsystems[1].level == 1;
        success *= top_data.subsystems[1].subsystems[0].id == "C";
        success *= top_data.subsystems[1].subsystems[0].level == 2;

        SystemT flat(readJson("SubsystemFlat.json").template get<DataT>());
        SystemT top(top_data);

        compareSystems(flat, top, success);

        return success.report(__func__);
      }

      /**
       * @brief Invalid hierarchies in the input data fail at construction.
       */
      TestOutcome invalidHierarchy()
      {
        TestStatus success = true;

        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::EVERYTHING);
        Log::misc() << "Testing that invalid hierarchies are rejected. "
                    << "Logged errors are expected.\n";
        Log::setVerbosity(previous_verbosity);

        const json hierarchical = readJson("SubsystemHierarchical.json");
        json       input;

        // Device uses a bus of a sibling subsystem
        input                                                  = hierarchical;
        input["subsystems"][0]["devices"][0]["ports"]["bus2"]  = 3;
        success                                               *= constructionThrows(input);

        // Device uses a signal of its parent
        const json signal_r = {{"signal_id", 1}, {"name", "sr"}};
        const json signal_i = {{"signal_id", 2}, {"name", "si"}};
        const json source   = {{"class", "ConstantSignalSource"},
                               {"ports", {{"sr", 1}, {"si", 2}}},
                               {"id", "1"},
                               {"params", json::object()}};
        input               = hierarchical;
        input["signals"]    = json::array({signal_r, signal_i});
        input["subsystems"][0]["devices"].push_back(source);
        success *= constructionThrows(input);

        // Devices in a system with subsystems
        input                = hierarchical;
        json load            = input["subsystems"][1]["subsystems"][0]["devices"][2];
        load["ports"]["bus"] = 1;
        input["devices"].push_back(load);
        success *= constructionThrows(input);

        // Programmatic hierarchy: devices and subsystems in one system
        {
          SystemT with_device(readJson("SubsystemFlat.json").template get<DataT>());
          SystemT with_subsystem;
          SystemT subsystem;
          SystemT other;
          with_subsystem.addSubsystem(&subsystem);

          bool add_subsystem_throws = false;
          try
          {
            with_device.addSubsystem(&other);
          }
          catch (const std::exception&)
          {
            add_subsystem_throws = true;
          }
          success *= add_subsystem_throws;

          bool add_component_throws = false;
          try
          {
            with_subsystem.addComponent(with_device.getComponent(0));
          }
          catch (const std::exception&)
          {
            add_component_throws = true;
          }
          success *= add_component_throws;
        }

        // Subsystem without an ID
        input = hierarchical;
        input["subsystems"][0].erase("id");
        success *= constructionThrows(input);

        // Subsystem without a level
        input = hierarchical;
        input["subsystems"][0].erase("level");
        success *= constructionThrows(input);

        // Subsystem level not matching the hierarchy (grandchild at level 1)
        input                                             = hierarchical;
        input["subsystems"][1]["subsystems"][0]["level"]  = 1;
        success                                          *= constructionThrows(input);

        return success.report(__func__);
      }

      /**
       * @brief Subsystems must have their parent's system bases, which
       * setSystemBase() passes down the hierarchy.
       */
      TestOutcome systemBases()
      {
        TestStatus success = true;

        // Programmatically built systems do not own their subsystems.
        SystemT top;
        SystemT a;
        SystemT b;
        top.addSubsystem(&a);

        b.setSystemBase(50.0, 100.0e6);
        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::EVERYTHING);
        Log::misc() << "Testing that a subsystem with different bases is rejected. "
                    << "Logged errors are expected.\n";
        Log::setVerbosity(previous_verbosity);
        try
        {
          a.addSubsystem(&b);
          success *= false;
        }
        catch (const std::exception&)
        {
        }

        // A gets the new bases from top, so B now fits.
        top.setSystemBase(50.0, 100.0e6);
        try
        {
          a.addSubsystem(&b);
        }
        catch (const std::exception&)
        {
          success *= false;
        }

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /**
       * @brief Compare the hierarchical Jacobian with the flat one and with
       * DependencyTracking, and the subsystem Jacobians with the top-level
       * one, before and after a fault changes the Jacobian values.
       */
      TestOutcome jacobian()
      {
        using DependencyMap = DependencyTracking::Variable::DependencyMap;

        TestStatus success = true;

        SystemT flat(readJson("SubsystemFlat.json").template get<DataT>());
        SystemT top(readJson("SubsystemHierarchical.json").template get<DataT>());
        success *= top.numSubsystems() == 2 && top.getSubsystem(1)->numSubsystems() == 1;
        if (top.numSubsystems() != 2 || top.getSubsystem(1)->numSubsystems() != 1)
        {
          return success.report(__func__);
        }
        auto* a = top.getSubsystem(0);
        auto* b = top.getSubsystem(1);
        auto* c = b->getSubsystem(0);

        PhasorDynamics::SystemModel<DependencyTracking::Variable, IdxT> top_dt(
            readJson("SubsystemHierarchical.json").template get<DataT>());

        for (auto* system : {&flat, &top})
        {
          success *= system->allocate() == 0;
          success *= system->initialize() == 0;
        }
        success *= top_dt.allocate() == 0;
        success *= top_dt.initialize() == 0;

        const IdxT n = top.size();

        // Duplicate entries may be summed in a different order.
        const RealT tol = 10 * std::numeric_limits<RealT>::epsilon();

        for (const bool fault : {false, true})
        {
          if (fault)
          {
            success *= flat.getBus(4)->setFault(true, 0.0, 0.01) == 0;
            success *= top.getBus(4)->setFault(true, 0.0, 0.01) == 0;
            success *= top_dt.getBus(4)->setFault(true, 0.0, 0.01) == 0;
          }
          for (auto* system : {&flat, &top})
          {
            success *= system->evaluateResidual() == 0;
            success *= system->evaluateJacobian() == 0;
          }
          success *= top_dt.evaluateResidual() == 0;
          success *= top_dt.evaluateJacobian() == 0;
          if (flat.getCsrJacobian() == nullptr || top.getCsrJacobian() == nullptr
              || top_dt.getCsrJacobian() == nullptr)
          {
            success *= false;
            return success.report(__func__);
          }

          // Both systems have the same ordering, see compareSystems().
          const auto jacobian_flat = MapFromCsr(flat.getCsrJacobian());
          const auto jacobian_top  = MapFromCsr(top.getCsrJacobian());
          const auto jacobian_dt   = MapFromCsr(top_dt.getCsrJacobian());
          for (IdxT i = 0; i < n; ++i)
          {
            success *= isEqual(jacobian_top[i], jacobian_flat[i]);
            success *= isEqual(jacobian_top[i], jacobian_dt[i], tol);
          }

          for (auto* subsystem : {a, b, c})
          {
            success *= subsystem->evaluateJacobian() == 0;
            if (subsystem->getCsrJacobian() == nullptr)
            {
              success *= false;
              return success.report(__func__);
            }

            // Subsystem Jacobian in global indices. Coupling with external
            // variables is not checked yet.
            const auto jacobian  = MapFromCsr(subsystem->getCsrJacobian());
            const IdxT begin     = subsystem->offset();
            const IdxT end       = begin + subsystem->size();
            success             *= static_cast<IdxT>(jacobian.size()) == end;
            if (static_cast<IdxT>(jacobian.size()) != end)
            {
              return success.report(__func__);
            }

            for (IdxT i = begin; i < end; ++i)
            {
              DependencyMap expected_row;
              DependencyMap actual_row;
              for (const auto& [column, value] : jacobian_top[i])
              {
                if (begin <= column && column < end)
                {
                  expected_row[column] = value;
                }
              }
              for (const auto& [column, value] : jacobian[i])
              {
                if (begin <= column && column < end)
                {
                  actual_row[column] = value;
                }
              }
              success *= isEqual(actual_row, expected_row, tol);
            }
          }
        }

        return success.report(__func__);
      }
#endif

    private:
      /// Read a JSON file copied next to the test executable
      static json readJson(const std::string& file_name)
      {
        std::ifstream file(file_name);
        return json::parse(file);
      }

      /// Whether constructing a system from the input throws
      static bool constructionThrows(const json& input)
      {
        try
        {
          const DataT data = input.template get<DataT>();
          SystemT     system(data);
        }
        catch (const std::exception&)
        {
          return true;
        }
        return false;
      }

      /**
       * @brief Compare top = [bus 1 | A[bus 2] | B[bus 3 | C[bus 4, load]]]
       * with the flat system, before and after a fault on bus 4.
       */
      static void compareSystems(SystemT& flat, SystemT& top, TestStatus& success)
      {
        for (auto* system : {&flat, &top})
        {
          success *= system->allocate() == 0;
          success *= system->initialize() == 0;
        }

        success *= top.numSubsystems() == 2;
        if (top.numSubsystems() != 2)
        {
          return;
        }
        auto* a  = top.getSubsystem(0);
        auto* b  = top.getSubsystem(1);
        success *= b->numSubsystems() == 1;
        if (b->numSubsystems() != 1)
        {
          return;
        }
        auto* c = b->getSubsystem(0);

        // Branches have no variables, so both systems order variables as
        // buses 1-4, then the load.
        success *= a->offset() == top.getBus(1)->size();
        success *= a->size() == a->getBus(2)->size();
        success *= b->offset() == a->offset() + a->size();
        success *= b->size() == b->getBus(3)->size() + c->size();
        success *= c->offset() == b->offset() + b->getBus(3)->size();
        success *= c->size() == c->getBus(4)->size() + c->getComponent(2)->size();
        success *= top.size() == top.getBus(1)->size() + a->size() + b->size();
        success *= top.size() == flat.size();
        if (top.size() != flat.size())
        {
          return;
        }

        // External buses
        BusT* bus1  = top.getBus(1);
        BusT* bus3  = b->getBus(3);
        success    *= a->getExternalBuses() == std::vector<BusT*>{bus1};
        success    *= b->getExternalBuses() == std::vector<BusT*>{bus1};
        success    *= c->getExternalBuses() == std::vector<BusT*>{bus1, bus3};

        // The top system finds subsystem buses.
        success *= top.getBus(3) == bus3;
        success *= top.getBus(4) == c->getBus(4);

        for (const bool fault : {false, true})
        {
          if (fault)
          {
            success *= flat.getBus(4)->setFault(true, 0.0, 0.01) == 0;
            success *= top.getBus(4)->setFault(true, 0.0, 0.01) == 0;
          }

          for (auto* system : {&flat, &top})
          {
            success *= system->evaluateResidual() == 0;
          }
          const auto* f_flat = flat.getResidual().getData();
          const auto* f_top  = top.getResidual().getData();
          for (IdxT i = 0; i < flat.size(); ++i)
          {
            success *= isEqual(f_top[i], f_flat[i]);
          }
        }
      }
    };
  } // namespace Testing
} // namespace GridKit
