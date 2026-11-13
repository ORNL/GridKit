/**
 * @file TwoAreaHierarchicalTests.hpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 * @brief Compare the TwoArea case with the same case split into hierarchical
 * subsystems, through an index map built from bus IDs and component labels.
 */
#pragma once

#include <limits>
#include <string>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class TwoAreaHierarchicalTests
    {
    private:
      using RealT         = ScalarT;
      using SystemT       = PhasorDynamics::SystemModel<ScalarT, IdxT>;
      using BusT          = PhasorDynamics::BusBase<ScalarT, IdxT>;
      using DependencyMap = DependencyTracking::Variable::DependencyMap;

    public:
      TwoAreaHierarchicalTests()  = default;
      ~TwoAreaHierarchicalTests() = default;

      /**
       * @brief Compare the flat and hierarchical TwoArea systems at the
       * initial point and after a short run with a fault.
       */
      TestOutcome flatVsHierarchical()
      {
        using namespace GridKit::PhasorDynamics;
        using namespace AnalysisManager::Sundials;
        using Log = GridKit::Utilities::Logger;

        TestStatus success = true;

        auto flat_data         = parseSystemModelData(std::string("TwoArea.case.json"));
        auto hierarchical_data = parseSystemModelData(std::string("TwoAreaHierarchical.case.json"));
        applyState("TwoArea.state.json", flat_data);
        applyState("TwoArea.state.json", hierarchical_data);

        SystemT flat(flat_data);
        SystemT hierarchical(hierarchical_data);

        success *= hierarchical.numSubsystems() == 2;
        if (hierarchical.numSubsystems() == 2)
        {
          auto* area1 = hierarchical.getSubsystem(0);
          auto* area2 = hierarchical.getSubsystem(1);

          // Areas have no subsystems
          success *= area1->numSubsystems() == 0;
          success *= area2->numSubsystems() == 0;

          for (auto* system : {&flat, &hierarchical})
          {
            success *= system->allocate() == 0;
            success *= system->initialize() == 0;
            success *= system->evaluateResidual() == 0;
            success *= system->evaluateJacobian() == 0;
          }

          const IdxT n  = flat.size();
          success      *= hierarchical.size() == n;
          if (hierarchical.size() == n)
          {
            // Hierarchical variable and residual index of each flat index
            std::vector<IdxT> variable_index(n, 0);
            std::vector<IdxT> residual_index(n, 0);
            for (IdxT bus_index = 0; bus_index < flat.numBuses(); ++bus_index)
            {
              const IdxT bus_id = flat.getBusId(bus_index);

              BusT* flat_bus = flat.getBus(bus_id);
              BusT* bus      = hierarchical.getBus(bus_id);
              for (IdxT variable = 0; variable < flat_bus->size(); ++variable)
              {
                const IdxT flat_variable      = flat_bus->getVariableIndex(variable);
                const IdxT flat_residual      = flat_bus->getResidualIndex(variable);
                variable_index[flat_variable] = bus->getVariableIndex(variable);
                residual_index[flat_residual] = bus->getResidualIndex(variable);
              }
            }
            for (IdxT component_index = 0; component_index < flat.numComponents(); ++component_index)
            {
              auto* flat_component = flat.getComponent(component_index);
              auto* component      = hierarchical.getComponent(flat.getComponentLabel(component_index));
              for (IdxT variable = 0; variable < flat_component->size(); ++variable)
              {
                const IdxT flat_variable      = flat_component->getVariableIndex(variable);
                const IdxT flat_residual      = flat_component->getResidualIndex(variable);
                variable_index[flat_variable] = component->getVariableIndex(variable);
                residual_index[flat_residual] = component->getResidualIndex(variable);
              }
            }

            // The map reaches every hierarchical index exactly once
            std::vector<IdxT> variable_count(n, 0);
            std::vector<IdxT> residual_count(n, 0);
            for (IdxT position = 0; position < n; ++position)
            {
              const IdxT variable = variable_index[position];
              const IdxT residual = residual_index[position];
              ++variable_count[variable];
              ++residual_count[residual];
            }
            for (IdxT position = 0; position < n; ++position)
            {
              success *= variable_count[position] == 1;
              success *= residual_count[position] == 1;
            }

            // External buses sum device contributions in a different order.
            const RealT residual_tol = 100 * std::numeric_limits<RealT>::epsilon();

            const auto* y_flat = flat.y().getData();
            const auto* y      = hierarchical.y().getData();
            const auto* f_flat = flat.getResidual().getData();
            const auto* f      = hierarchical.getResidual().getData();
            for (IdxT position = 0; position < n; ++position)
            {
              const IdxT variable  = variable_index[position];
              const IdxT residual  = residual_index[position];
              success             *= isEqual(y[variable], y_flat[position]);
              success             *= isEqual(f[residual], f_flat[position], residual_tol);
            }

#ifdef GRIDKIT_ENABLE_ENZYME
            auto* jacobian_flat  = flat.getCsrJacobian();
            auto* jacobian       = hierarchical.getCsrJacobian();
            success             *= jacobian_flat != nullptr && jacobian != nullptr;
            if (jacobian_flat != nullptr && jacobian != nullptr)
            {
              success *= jacobian->getNnz() == jacobian_flat->getNnz();

              const auto expected = MapFromCsr(jacobian);

              // Flat Jacobian entries at their hierarchical row and column
              const auto                 flat_entries = MapFromCsr(jacobian_flat);
              std::vector<DependencyMap> mapped_flat(n);
              for (IdxT row = 0; row < n; ++row)
              {
                DependencyMap& mapped_row = mapped_flat[residual_index[row]];
                for (const auto& [column, value] : flat_entries[row])
                {
                  const IdxT mapped_column  = variable_index[column];
                  mapped_row[mapped_column] = value;
                }
              }
              for (IdxT row = 0; row < n; ++row)
              {
                success *= isEqual(expected[row], mapped_flat[row]);
              }

              // Bus 8 entries plus the Jacobians of both areas give the hierarchical
              // system Jacobian.
              BusT*                      bus8 = hierarchical.getBus(8);
              std::vector<DependencyMap> sum(n);
              auto*                      bus_jacobian = bus8->getCooJacobian();
              for (IdxT entry = 0; entry < bus_jacobian->getNnz(); ++entry)
              {
                const IdxT     row      = bus_jacobian->getRowData()[entry];
                const IdxT     column   = bus_jacobian->getColData()[entry];
                DependencyMap& sum_row  = sum[row];
                sum_row[column]        += bus_jacobian->getValues()[entry];
              }

              for (auto* area : {area1, area2})
              {
                success *= area->evaluateJacobian() == 0;
                success *= area->getExternalBuses() == std::vector<BusT*>{bus8};
                success *= area->getCsrJacobian() != nullptr;
                if (area->getCsrJacobian() != nullptr)
                {
                  const auto area_jacobian = MapFromCsr(area->getCsrJacobian());
                  for (IdxT row = 0; row < area_jacobian.size(); ++row)
                  {
                    DependencyMap& sum_row = sum[row];
                    for (const auto& [column, value] : area_jacobian[row])
                    {
                      sum_row[column] += value;
                    }
                  }
                }
              }

              // Duplicate entries may be summed in a different order.
              const RealT jacobian_tol = 10 * std::numeric_limits<RealT>::epsilon();
              for (IdxT row = 0; row < n; ++row)
              {
                success *= isEqual(sum[row], expected[row], jacobian_tol);
              }
            }
#endif

            // Short run with a fault at bus 7 in Area 1, set up as in the
            // DynamicSimulation application with the example studies' settings
            const RealT dt_monitor   = 1.0 / 240.0;
            const RealT relative_tol = 1.0e-7;
            const RealT absolute_tol = 1.0e-9;
            const RealT fault_x      = 0.01;
            for (auto* system : {&flat, &hierarchical})
            {
              const char* name = "hierarchical";
              if (system == &flat)
              {
                name = "flat";
              }
              BusT* bus7 = system->getBus(7);

              try
              {
                Ida<ScalarT, IdxT> ida(system);
                ida.setTolerance(relative_tol, absolute_tol);
                ida.setConsistentICType(IdaConsistentICType::YA_YDP);
                ida.configureSimulation();
                ida.initializeSimulation(0.0);
                ida.runSimulation(0.1, dt_monitor);

                success *= bus7->setFault(true, 0.0, fault_x) == 0;
                ida.initializeSimulation(0.1);
                ida.runSimulation(0.2, dt_monitor);

                success *= bus7->setFault(false, 0.0, fault_x) == 0;
                ida.initializeSimulation(0.2);
                ida.runSimulation(1.0, dt_monitor);
              }
              catch (const SundialsException&)
              {
                Log::error() << "IDA failed for the " << name << " system\n";
                success *= false;
              }
            }

            // The two systems order variables and sum residuals differently, so
            // IDA takes different steps. Voltages and states agree to a multiple
            // of the integrator tolerance, not to eps.
            const RealT voltage_tol = 10 * relative_tol;
            for (IdxT bus_index = 0; bus_index < flat.numBuses(); ++bus_index)
            {
              const IdxT  bus_id    = flat.getBusId(bus_index);
              const BusT* flat_bus  = flat.getBus(bus_id);
              const BusT* bus       = hierarchical.getBus(bus_id);
              success              *= isEqual(bus->Vr(), flat_bus->Vr(), voltage_tol);
              success              *= isEqual(bus->Vi(), flat_bus->Vi(), voltage_tol);
            }

            const RealT state_tol = 100 * relative_tol;
            for (IdxT position = 0; position < n; ++position)
            {
              const IdxT variable  = variable_index[position];
              success             *= isEqual(y[variable], y_flat[position], state_tol);
            }
          }
        }

        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
