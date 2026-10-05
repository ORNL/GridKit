#pragma once

#include <cstddef>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <GridKit/Model/Coupling.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// One entry of a .partition.json file
    struct PartitionData
    {
      std::string              name;
      std::vector<std::size_t> buses;
    };

    /// A case split into one case per partition, with maps back to the case
    struct PartitionedModelData
    {
      std::vector<SystemModelData<>>                   partitions;    ///< One case per partition
      std::map<std::size_t, std::size_t>               bus_partition; ///< Bus ID -> owning partition
      std::vector<std::pair<std::size_t, std::size_t>> components;    ///< Original ID -> (partition, ID)
      std::vector<std::pair<std::size_t, std::size_t>> faults;        ///< Original ID -> (partition, ID)
    };

    std::vector<PartitionData> parsePartitionData(const std::filesystem::path& file);

    /**
     * @brief Split a case into one case per partition.
     *
     * Buses and devices go to the partition owning their bus; busless
     * controllers follow their signals. A branch between partitions is copied
     * to both sides, and its far end becomes a BusInfinite external bus there.
     */
    PartitionedModelData partitionSystemModelData(const SystemModelData<>&          model,
                                                  const std::vector<PartitionData>& partitions);

    /**
     * @brief Couple each partition's external buses to the partition that
     * solves for them.
     *
     * An external bus no partition solves for (an actual infinite bus) keeps
     * its prescribed voltage.
     */
    template <class ScalarT, typename IdxT>
    std::vector<std::vector<Model::Coupling<ScalarT, IdxT>>>
    connectPartitions(const std::vector<std::unique_ptr<SystemModel<ScalarT, IdxT>>>& partitions)
    {
      std::vector<std::vector<Model::Coupling<ScalarT, IdxT>>> couplings(partitions.size());
      for (std::size_t p = 0; p < partitions.size(); ++p)
      {
        for (auto* external : partitions[p]->externalBuses())
        {
          for (const auto& source : partitions)
          {
            auto* bus = source->findBus(external->busID());
            if (bus != nullptr && bus->size() > 0)
            {
              couplings[p].push_back({source.get(), bus->getVariableIndices()[0], &external->VrInput()});
              couplings[p].push_back({source.get(), bus->getVariableIndices()[1], &external->ViInput()});
            }
          }
        }
      }
      return couplings;
    }
  } // namespace PhasorDynamics
} // namespace GridKit
