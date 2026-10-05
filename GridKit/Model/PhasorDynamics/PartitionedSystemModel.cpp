#include "PartitionedSystemModel.hpp"

#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace
    {
      using Log   = ::GridKit::Utilities::Logger;
      using DataT = SystemModelData<double, size_t>;
      using IdxT  = size_t;

      void require(bool valid, const std::string& message)
      {
        if (!valid)
        {
          Log::error() << "Partition: " << message << std::endl;
          throw std::runtime_error("Partition: " + message);
        }
      }

      template <typename Device, typename Visit>
      void forEachSignal(const Device& device, Visit&& visit)
      {
        for (const auto& [port, signal] : device.signal_inputs)
        {
          visit(signal);
        }
        for (const auto& [port, signal] : device.signal_outputs)
        {
          visit(signal);
        }
      }

      /// A device sits in the partition of its buses; a device without a bus
      /// (governor, stabilizer) in the partition of a signal it shares.
      template <typename Device>
      std::optional<IdxT> devicePartition(const Device&               device,
                                          const std::map<IdxT, IdxT>& bus_partition,
                                          const std::map<IdxT, IdxT>& signal_partition)
      {
        std::optional<IdxT> partition;
        for (const auto& [port, bus] : device.buses)
        {
          require(bus_partition.contains(bus), device.disambiguation_string + " connects to unknown bus " + std::to_string(bus));
          const auto p = bus_partition.at(bus);
          require(!partition || *partition == p, device.disambiguation_string + " connects two partitions");
          partition = p;
        }
        forEachSignal(device, [&](IdxT signal)
                      {
          if (!partition && signal_partition.contains(signal))
          {
            partition = signal_partition.at(signal);
          } });
        return partition;
      }

      /// Signals follow the devices that use them. Control loops are never split.
      std::map<IdxT, IdxT> signalPartitions(const DataT& model, const std::map<IdxT, IdxT>& bus_partition)
      {
        std::map<IdxT, IdxT> signal_partition;
        bool                 changed = true;
        while (changed)
        {
          changed = false;
          forEachDeviceGroup<DataT>([&](auto member)
                                    {
            using Device = typename std::decay_t<decltype(model.*member)>::value_type;
            if constexpr (!std::is_same_v<Device, DataT::BranchDataT>)
            {
              for (const auto& device : model.*member)
              {
                const auto p = devicePartition(device, bus_partition, signal_partition);
                if (!p)
                {
                  continue;
                }
                forEachSignal(device, [&](IdxT signal)
                              {
                  const auto [entry, inserted] = signal_partition.emplace(signal, *p);
                  require(entry->second == *p, "signal " + std::to_string(signal) + " is used in two partitions");
                  changed |= inserted; });
              }
            } });
        }
        return signal_partition;
      }

      /// One system per partition, including each tie branch's far end.
      std::vector<DataT> splitSystemModelData(const DataT&                                      model,
                                              const std::vector<PartitionData<double, size_t>>& partitions,
                                              const std::map<IdxT, IdxT>&                       bus_partition)
      {
        const auto n = partitions.size();
        require(n >= 2, "at least two partitions are required");
        std::vector<DataT> systems(n);
        for (IdxT p = 0; p < n; ++p)
        {
          auto& system            = systems[p];
          system.format_version   = model.format_version;
          system.format_revision  = model.format_revision;
          system.case_name        = model.case_name + ":" + partitions[p].name;
          system.case_date_time   = model.case_date_time;
          system.case_description = model.case_description;
          system.case_comments    = model.case_comments;
          system.freq_base        = model.freq_base;
          system.va_base          = model.va_base;
        }

        std::map<IdxT, const DataT::BusDataT*> case_bus;
        for (const auto& bus : model.bus)
        {
          require(bus_partition.contains(bus.bus_id), "bus " + std::to_string(bus.bus_id) + " is in no partition");
          systems[bus_partition.at(bus.bus_id)].bus.push_back(bus);
          case_bus.emplace(bus.bus_id, &bus);
        }
        require(case_bus.size() == bus_partition.size(), "partition file names a bus that is not in the case");

        const auto signal_partition = signalPartitions(model, bus_partition);
        for (const auto& signal : model.signal)
        {
          if (signal_partition.contains(signal.signal_id))
          {
            systems[signal_partition.at(signal.signal_id)].signal.push_back(signal);
          }
        }

        std::vector<std::set<IdxT>> far_ends(n);
        forEachDeviceGroup<DataT>([&](auto member)
                                  {
          using Device = typename std::decay_t<decltype(model.*member)>::value_type;
          for (const auto& device : model.*member)
          {
            if constexpr (std::is_same_v<Device, DataT::BranchDataT>)
            {
              // A tie branch is in both partitions; the bus1 side keeps the monitors.
              const auto a = device.buses.at(BranchBuses::bus1);
              const auto b = device.buses.at(BranchBuses::bus2);
              const auto p = bus_partition.at(a);
              const auto q = bus_partition.at(b);
              systems[p].branch.push_back(device);
              if (q != p)
              {
                auto copy = device;
                copy.monitored_variables.clear();
                systems[q].branch.push_back(copy);
                far_ends[p].insert(b);
                far_ends[q].insert(a);
              }
            }
            else
            {
              const auto p = devicePartition(device, bus_partition, signal_partition);
              require(p.has_value(), device.disambiguation_string + " has no bus or signal");
              (systems[*p].*member).push_back(device);
            }
          } });

        // Each tie branch's far end: the other partition's bus at its case voltage.
        for (IdxT p = 0; p < n; ++p)
        {
          for (const auto id : far_ends[p])
          {
            auto bus     = *case_bus.at(id);
            bus.bus_type = DataT::BusDataT::BusType::SLACK;
            bus.monitored_variables.clear();
            systems[p].bus.push_back(bus);
          }
        }
        return systems;
      }
    } // namespace

    template <typename scalar_type, typename index_type>
    PartitionedSystemModel<scalar_type, index_type>::PartitionedSystemModel(const SystemModelData<RealT, IdxT>&            data,
                                                                            const std::vector<PartitionData<RealT, IdxT>>& partitions)
    {
      for (IdxT p = 0; p < partitions.size(); ++p)
      {
        for (const auto bus : partitions[p].buses)
        {
          require(bus_partition_.emplace(bus, p).second, "bus " + std::to_string(bus) + " is in two partitions");
        }
      }
      for (const auto& system : splitSystemModelData(data, partitions, bus_partition_))
      {
        partitions_.push_back(std::make_unique<SystemModelT>(system));
      }
    }

    template <typename scalar_type, typename index_type>
    typename PartitionedSystemModel<scalar_type, index_type>::IdxT
    PartitionedSystemModel<scalar_type, index_type>::numPartitions() const
    {
      return static_cast<IdxT>(partitions_.size());
    }

    template <typename scalar_type, typename index_type>
    typename PartitionedSystemModel<scalar_type, index_type>::SystemModelT&
    PartitionedSystemModel<scalar_type, index_type>::getPartition(IdxT partition)
    {
      return *partitions_.at(partition);
    }

    template <typename scalar_type, typename index_type>
    typename PartitionedSystemModel<scalar_type, index_type>::BusT*
    PartitionedSystemModel<scalar_type, index_type>::getBus(IdxT bus_id)
    {
      return partitions_[bus_partition_.at(bus_id)]->getBus(bus_id);
    }

    template class PartitionedSystemModel<double, size_t>;
  } // namespace PhasorDynamics
} // namespace GridKit
