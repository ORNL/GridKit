/**
 * @file PartitionPlacement.hpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 * @brief Place the contents of flat system model data in the partitions of a partition file.
 */
#pragma once

#include <map>
#include <string>
#include <vector>

#include <GridKit/Model/ConfigurationChecks.hpp>
#include <GridKit/Model/PhasorDynamics/Partition/PartitionData.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Generates hierarchical system model data from flat data and a partition
     *
     * The system tree follows the partition tree. Each partition becomes a
     * subsystem, and each leaf gets its internal buses and the devices placed
     * in it. External buses go to the top system. Case metadata and monitor
     * sinks stay in the top system. The placement rules are in README.md in
     * this directory.
     */
    template <typename real_type, typename index_type>
    class PartitionPlacement
    {
      using RealT          = real_type;
      using IdxT           = index_type;
      using DataT          = SystemModelData<RealT, IdxT>;
      using PartitionDataT = PartitionData<IdxT>;

    public:
      /**
       * @brief Build the system tree of `partition` from flat `data`
       *
       * @throws std::runtime_error for a badly formed partition
       */
      static DataT apply(const DataT& data, const PartitionDataT& partition);

    private:
      PartitionPlacement(const DataT& data);

      DataT place(const PartitionDataT& partition);
      void  buildSystem(const PartitionDataT& partition, DataT& system);
      void  addBusList(const std::vector<IdxT>& buses,
                       const std::string&       name,
                       bool                     internal,
                       IdxT                     leaf);

      const DataT&               data_;   ///< Flat input data
      Model::ConfigurationChecks checks_; ///< Errors of a badly formed partition

      std::map<IdxT, IdxT> bus_index_;    ///< Bus ID -> index in data_.bus
      std::map<IdxT, IdxT> signal_bus_;   ///< Signal ID -> index of its bus
      std::map<IdxT, IdxT> signal_group_; ///< Signal ID -> group of the devices using it

      std::vector<std::string>       device_labels_; ///< Label of each device, for messages
      std::vector<std::vector<IdxT>> device_buses_;  ///< Bus indices of each device
      std::vector<IdxT>              device_group_;  ///< Group of each device, its lowest index
      std::vector<IdxT>              next_in_group_; ///< Next device of its group, or INVALID_INDEX
      std::vector<IdxT>              group_leaf_;    ///< Leaf of each group, at its first device

      std::vector<IdxT> bus_internal_leaf_; ///< Leaf each bus is internal to, or INVALID_INDEX
      std::vector<IdxT> bus_last_leaf_;     ///< Last leaf listing each bus, or INVALID_INDEX
      IdxT              num_leaves_{0};     ///< Leaves numbered so far, in file order
    };
  } // namespace PhasorDynamics
} // namespace GridKit
