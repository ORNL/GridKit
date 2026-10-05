#pragma once

#include <map>
#include <memory>
#include <vector>

#include <GridKit/Model/PhasorDynamics/PartitionData.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief A system model split into one system model per partition.
     *
     * Each partition solves for its own buses and holds every device on them;
     * devices without a bus follow their signals. A tie branch is in both
     * partitions it joins; its far end is the other partition's bus, held at
     * its case voltage.
     */
    template <typename scalar_type, typename index_type>
    class PartitionedSystemModel
    {
    public:
      using ScalarT      = scalar_type;
      using IdxT         = index_type;
      using SystemModelT = SystemModel<ScalarT, IdxT>;
      using RealT        = typename SystemModelT::RealT;
      using BusT         = typename SystemModelT::BusT;

      PartitionedSystemModel(const SystemModelData<RealT, IdxT>&            data,
                             const std::vector<PartitionData<RealT, IdxT>>& partitions);

      IdxT          numPartitions() const;
      SystemModelT& getPartition(IdxT partition);

      /// The bus with this ID in the partition that solves for it
      BusT* getBus(IdxT bus_id);

    private:
      std::vector<std::unique_ptr<SystemModelT>> partitions_;
      std::map<IdxT, IdxT>                       bus_partition_; ///< Partition that solves for each bus, by bus ID
    };
  } // namespace PhasorDynamics
} // namespace GridKit
