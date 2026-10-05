#include "PartitionData.hpp"

#include <fstream>
#include <optional>
#include <set>
#include <stdexcept>
#include <type_traits>

#include <nlohmann/json.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace
    {
      using DataT = SystemModelData<>;
      using IdxT  = std::size_t;

      void require(bool valid, const std::string& message)
      {
        if (!valid)
        {
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

      /// A device sits in the partition of its buses; a busless device
      /// (governor, stabilizer) in the partition of a signal it shares.
      template <typename Device>
      std::optional<std::size_t> devicePartition(const Device&                      device,
                                                 const std::map<IdxT, std::size_t>& bus_partition,
                                                 const std::map<IdxT, std::size_t>& signal_partition)
      {
        std::optional<std::size_t> partition;
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
      std::map<IdxT, std::size_t> signalPartitions(const DataT& model, const std::map<IdxT, std::size_t>& bus_partition)
      {
        std::map<IdxT, std::size_t> signal_partition;
        bool                        changed = true;
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
    } // namespace

    std::vector<PartitionData> parsePartitionData(const std::filesystem::path& file)
    {
      std::ifstream stream(file);
      require(stream.good(), "could not open " + file.string());
      const auto document = nlohmann::json::parse(stream);

      std::vector<PartitionData> partitions;
      for (const auto& entry : document.at("partitions"))
      {
        partitions.push_back({entry.at("name").get<std::string>(),
                              entry.at("buses").get<std::vector<std::size_t>>()});
      }
      return partitions;
    }

    PartitionedModelData partitionSystemModelData(const DataT& model, const std::vector<PartitionData>& partitions)
    {
      const auto           n = partitions.size();
      PartitionedModelData result;
      require(n >= 2, "at least two partitions are required");
      result.partitions.resize(n);

      for (std::size_t p = 0; p < n; ++p)
      {
        auto& data            = result.partitions[p];
        data.format_version   = model.format_version;
        data.format_revision  = model.format_revision;
        data.case_name        = model.case_name + ":" + partitions[p].name;
        data.case_date_time   = model.case_date_time;
        data.case_description = model.case_description;
        data.case_comments    = model.case_comments;
        data.freq_base        = model.freq_base;
        data.va_base          = model.va_base;
        for (const auto bus : partitions[p].buses)
        {
          require(result.bus_partition.emplace(bus, p).second, "bus " + std::to_string(bus) + " is in two partitions");
        }
      }

      std::map<IdxT, const DataT::BusDataT*> case_bus;
      for (const auto& bus : model.bus)
      {
        require(result.bus_partition.contains(bus.bus_id), "bus " + std::to_string(bus.bus_id) + " is in no partition");
        result.partitions[result.bus_partition.at(bus.bus_id)].bus.push_back(bus);
        case_bus.emplace(bus.bus_id, &bus);
      }
      require(case_bus.size() == result.bus_partition.size(), "partition file names a bus that is not in the case");

      const auto signal_partition = signalPartitions(model, result.bus_partition);
      for (const auto& signal : model.signal)
      {
        if (signal_partition.contains(signal.signal_id))
        {
          result.partitions[signal_partition.at(signal.signal_id)].signal.push_back(signal);
        }
      }

      std::vector<IdxT>           next_component(n, 0);
      std::vector<IdxT>           next_fault(n, 0);
      std::vector<std::set<IdxT>> external(n);
      forEachDeviceGroup<DataT>([&](auto member)
                                {
        using Device = typename std::decay_t<decltype(model.*member)>::value_type;
        for (const auto& device : model.*member)
        {
          if constexpr (std::is_same_v<Device, DataT::BranchDataT>)
          {
            // Each side evaluates its own copy and keeps only its own end; the
            // far end is an external bus. The bus1 side keeps the monitors.
            const auto a = device.buses.at(BranchBuses::bus1);
            const auto b = device.buses.at(BranchBuses::bus2);
            const auto p = result.bus_partition.at(a);
            const auto q = result.bus_partition.at(b);
            result.partitions[p].branch.push_back(device);
            result.components.emplace_back(p, next_component[p]++);
            if (q != p)
            {
              auto copy = device;
              copy.monitored_variables.clear();
              result.partitions[q].branch.push_back(copy);
              ++next_component[q];
              external[p].insert(b);
              external[q].insert(a);
            }
          }
          else
          {
            const auto p = devicePartition(device, result.bus_partition, signal_partition);
            require(p.has_value(), device.disambiguation_string + " has no bus or signal");
            (result.partitions[*p].*member).push_back(device);
            result.components.emplace_back(*p, next_component[*p]++);
            if constexpr (std::is_same_v<Device, DataT::BusFaultDataT>)
            {
              result.faults.emplace_back(*p, next_fault[*p]++);
            }
          }
        } });

      // An external bus is another partition's bus whose voltage is an input here.
      for (std::size_t p = 0; p < n; ++p)
      {
        for (const auto id : external[p])
        {
          auto bus     = *case_bus.at(id);
          bus.bus_type = DataT::BusDataT::BusType::SLACK;
          bus.monitored_variables.clear();
          result.partitions[p].bus.push_back(bus);
        }
      }
      return result;
    }
  } // namespace PhasorDynamics
} // namespace GridKit
