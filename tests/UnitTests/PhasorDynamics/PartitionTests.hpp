#pragma once

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include <GridKit/Model/PhasorDynamics/BusFault/BusFault.hpp>
#include <GridKit/Model/PhasorDynamics/PartitionData.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class PartitionTests
    {
      using SystemT = PhasorDynamics::SystemModel<ScalarT, IdxT>;
      using RealT   = typename SystemT::RealT;

    public:
      /// The partitions reproduce the intact case's residual at a perturbed
      /// state with every fault applied.
      TestOutcome residualMatchesIntactCase(const std::string& case_file, const std::string& partition_file)
      {
        TestStatus success = true;

        const auto data  = PhasorDynamics::parseSystemModelData(case_file);
        const auto parts = PhasorDynamics::partitionSystemModelData(data, PhasorDynamics::parsePartitionData(partition_file));

        SystemT intact(data);
        intact.allocate();
        perturb(intact.y());
        perturb(intact.yp());

        std::vector<std::unique_ptr<SystemT>> systems;
        for (const auto& partition : parts.partitions)
        {
          systems.emplace_back(std::make_unique<SystemT>(partition))->allocate();
        }
        for (const auto& bus : data.bus)
        {
          copyState(*intact.getBus(bus.bus_id), *systems[parts.bus_partition.at(bus.bus_id)]->getBus(bus.bus_id));
        }
        for (IdxT id = 0; id < parts.components.size(); ++id)
        {
          const auto& [p, local] = parts.components[id];
          copyState(*intact.getComponent(id), *systems[p]->getComponent(local));
        }
        for (IdxT fault = 0; fault < parts.faults.size(); ++fault)
        {
          const auto& [p, local] = parts.faults[fault];
          intact.getBusFault(fault)->setStatus(true);
          systems[p]->getBusFault(local)->setStatus(true);
        }
        for (const auto& couplings : PhasorDynamics::connectPartitions(systems))
        {
          for (const auto& coupling : couplings)
          {
            coupling.input->value = coupling.source->y().getData()[coupling.index];
          }
        }

        intact.updateTime(0.0, 0.0);
        intact.evaluateResidual();
        for (auto& system : systems)
        {
          system->updateTime(0.0, 0.0);
          system->evaluateResidual();
        }

        for (const auto& bus : data.bus)
        {
          success *= sameResidual(*intact.getBus(bus.bus_id), *systems[parts.bus_partition.at(bus.bus_id)]->getBus(bus.bus_id));
        }
        for (IdxT id = 0; id < parts.components.size(); ++id)
        {
          const auto& [p, local]  = parts.components[id];
          success                *= sameResidual(*intact.getComponent(id), *systems[p]->getComponent(local));
        }

        return success.report(__func__);
      }

    private:
      /// Deterministic perturbation that also moves zero entries.
      static void perturb(typename SystemT::VectorT& v)
      {
        auto* data = v.getData();
        for (IdxT i = 0; i < v.getSize(); ++i)
        {
          data[i] = data[i] * (1.0 + 1.0e-2 * std::sin(i + 1.0)) + 1.0e-3 * std::cos(i + 1.0);
        }
        v.setDataUpdated();
      }

      template <class ModelT>
      static void copyState(ModelT& source, ModelT& target)
      {
        std::copy_n(source.y().getData(), source.size(), target.y().getData());
        std::copy_n(source.yp().getData(), source.size(), target.yp().getData());
        target.y().setDataUpdated();
        target.yp().setDataUpdated();
      }

      template <class ModelT>
      static bool sameResidual(ModelT& a, ModelT& b)
      {
        bool        same = true;
        const auto* fa   = a.getResidual().getData();
        const auto* fb   = b.getResidual().getData();
        for (IdxT i = 0; i < a.size(); ++i)
        {
          if (!isEqual(fb[i], fa[i]))
          {
            same = false;
          }
        }
        return same;
      }
    };
  } // namespace Testing
} // namespace GridKit
