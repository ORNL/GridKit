/**
 * @file PartitionPlacement.cpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 * @brief Place the contents of flat system model data in the partitions of a partition file.
 */
#include "PartitionPlacement.hpp"

#include <format>
#include <stdexcept>

#include <GridKit/Constants.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using Log = ::GridKit::Utilities::Logger;

    /**
     * @brief Empty one device vector of a system
     *
     * `forEachDeviceList` calls this once per device type. The top system
     * starts as a copy of the flat case, so its devices are cleared.
     */
    struct ClearDevices
    {
      template <class DevicesT>
      void operator()(DevicesT& devices)
      {
        devices.clear();
      }
    };

    /**
     * @brief Record the label, buses and signals of each device
     *
     * `forEachDeviceList` calls this once per device type. Each device appends
     * one entry to `device_labels` and `device_buses`, so devices are numbered
     * in `forEachDeviceList` order. Buses are stored as case bus indices.
     *
     * Each signal port of a device appends one entry to `signal_use_device`
     * (the device index) and `signal_use_id` (the signal ID).
     *
     * @todo Place devices using a bus signal in the leaf where the bus is internal.
     */
    template <typename index_type>
    struct CollectDevices
    {
      using IdxT = index_type;

      const std::map<IdxT, IdxT>&     bus_index;
      const std::map<IdxT, IdxT>&     signal_bus;
      std::vector<std::string>&       device_labels;
      std::vector<std::vector<IdxT>>& device_buses;
      std::vector<IdxT>&              signal_use_device;
      std::vector<IdxT>&              signal_use_id;
      Model::ConfigurationChecks&     checks;

      template <class DeviceDataT>
      void operator()(const std::vector<DeviceDataT>& devices)
      {
        for (const auto& device_data : devices)
        {
          const IdxT device = device_labels.size();
          device_labels.push_back(device_data.label());
          device_buses.emplace_back();

          // Case bus indices of the device's buses
          for (const auto& port : device_data.buses)
          {
            checks.check(bus_index.contains(port.second),
                         std::format("device {} uses unknown bus {}",
                                     device_data.label(),
                                     port.second));
            if (bus_index.contains(port.second))
            {
              device_buses[device].push_back(bus_index.at(port.second));
            }
          }

          // A signal of a bus stays in the bus's system, and subsystems do not
          // see their parent's signals. Devices using one are rejected for now.
          for (const auto& port : device_data.signal_inputs)
          {
            signal_use_device.push_back(device);
            signal_use_id.push_back(port.second);
            checks.check(!signal_bus.contains(port.second),
                         std::format("device {} uses signal {} of a bus",
                                     device_data.label(),
                                     port.second));
          }
          for (const auto& port : device_data.signal_outputs)
          {
            signal_use_device.push_back(device);
            signal_use_id.push_back(port.second);
            checks.check(!signal_bus.contains(port.second),
                         std::format("device {} uses signal {} of a bus",
                                     device_data.label(),
                                     port.second));
          }
        }
      }
    };

    /**
     * @brief Copy the devices placed in one leaf from the flat case
     *
     * `forEachDeviceList` calls this once per device type. It passes the flat
     * case's vector and the leaf system's vector of the same type.
     *
     * `device` numbers the devices of all types in `forEachDeviceList` order, as
     * `device_group` does. It starts at 0 for each leaf.
     */
    template <typename index_type>
    struct CopyLeafDevices
    {
      using IdxT = index_type;

      const std::vector<IdxT>& device_group;
      const std::vector<IdxT>& group_leaf;
      IdxT                     leaf;
      IdxT                     device;

      template <class DeviceDataT>
      void operator()(const std::vector<DeviceDataT>& devices,
                      std::vector<DeviceDataT>&       system_devices)
      {
        for (const auto& device_data : devices)
        {
          const IdxT group = device_group[device];
          if (group_leaf[group] == leaf)
          {
            system_devices.push_back(device_data);
          }
          ++device;
        }
      }
    };

    template <typename real_type, typename index_type>
    SystemModelData<real_type, index_type> PartitionPlacement<real_type, index_type>::apply(
        const DataT&          data,
        const PartitionDataT& partition)
    {
      if (!data.subsystems.empty())
      {
        Log::error() << "A partition file can only be applied to a flat case\n";
        throw std::runtime_error("PartitionPlacement::apply failed");
      }

      PartitionPlacement placement(data);
      return placement.place(partition);
    }

    /**
     * Records the buses and signal groups of the devices. They do not depend
     * on the partition.
     */
    template <typename real_type, typename index_type>
    PartitionPlacement<real_type, index_type>::PartitionPlacement(const DataT& data)
      : data_(data),
        checks_("Partition file")
    {
      // Buses and their signals
      for (IdxT bus = 0; bus < data_.bus.size(); ++bus)
      {
        const auto& bus_data = data_.bus[bus];
        bus_index_.emplace(bus_data.bus_id, bus);
        for (const auto& port : bus_data.signal_inputs)
        {
          signal_bus_[port.second] = bus;
        }
        for (const auto& port : bus_data.signal_outputs)
        {
          signal_bus_[port.second] = bus;
        }
      }
      bus_internal_leaf_.assign(data_.bus.size(), INVALID_INDEX<IdxT>);
      bus_last_leaf_.assign(data_.bus.size(), INVALID_INDEX<IdxT>);

      std::vector<IdxT>    signal_use_device;
      std::vector<IdxT>    signal_use_id;
      CollectDevices<IdxT> collect_devices{bus_index_,
                                           signal_bus_,
                                           device_labels_,
                                           device_buses_,
                                           signal_use_device,
                                           signal_use_id,
                                           checks_};
      DataT::forEachDeviceList(collect_devices, data_);

      // Label each device with the lowest device index of its signal group.
      const IdxT num_devices = device_labels_.size();
      device_group_.resize(num_devices);

      // Each device starts as its own group
      for (IdxT device = 0; device < num_devices; ++device)
      {
        device_group_[device] = device;
      }

      // Signals chain devices together, so labels are passed along until none
      // changes.
      bool changed = true;
      while (changed)
      {
        changed = false;
        for (IdxT use = 0; use < signal_use_id.size(); ++use)
        {
          const IdxT device    = signal_use_device[use];
          const IdxT signal_id = signal_use_id[use];
          signal_group_.emplace(signal_id, device_group_[device]);
          const IdxT group = signal_group_.at(signal_id);
          if (group < device_group_[device])
          {
            device_group_[device] = group;
            changed               = true;
          }
          else if (device_group_[device] < group)
          {
            signal_group_[signal_id] = device_group_[device];
            changed                  = true;
          }
        }
      }

      // Devices of each signal group, chained from its first device
      next_in_group_.assign(num_devices, INVALID_INDEX<IdxT>);
      // Last device chained so far, at the first device of each group
      std::vector<IdxT> last_in_group(num_devices);

      // Chain each device after the last device of its group
      for (IdxT device = 0; device < num_devices; ++device)
      {
        const IdxT group = device_group_[device];
        if (group != device)
        {
          next_in_group_[last_in_group[group]] = device;
        }
        last_in_group[group] = device;
      }

      group_leaf_.assign(num_devices, INVALID_INDEX<IdxT>);
    }

    /**
     * @brief Build the system tree of `partition` and check the placement
     *
     * Every case bus must be listed and every signal group placed in a leaf.
     * External buses and the signals of no leaf stay in the top system.
     *
     * @throws std::runtime_error for a badly formed partition
     */
    template <typename real_type, typename index_type>
    SystemModelData<real_type, index_type> PartitionPlacement<real_type, index_type>::place(
        const PartitionDataT& partition)
    {
      // Case metadata and monitor sinks stay in the top system.
      DataT top = data_;
      top.bus.clear();
      top.signal.clear();
      ClearDevices clear_devices;
      DataT::forEachDeviceList(clear_devices, top);

      buildSystem(partition, top);

      // Every case bus must be listed. A repeated case bus is never listed, so
      // it is reported here too.
      for (IdxT bus = 0; bus < data_.bus.size(); ++bus)
      {
        checks_.check(bus_last_leaf_[bus] != INVALID_INDEX<IdxT>,
                      std::format("bus {} is missing", data_.bus[bus].bus_id));
      }

      // Every signal group must be placed in a leaf
      for (IdxT device = 0; device < device_group_.size(); ++device)
      {
        if (device_group_[device] == device)
        {
          checks_.check(group_leaf_[device] != INVALID_INDEX<IdxT>,
                        std::format("no partition holds all buses of device {} and the devices"
                                    " connected to it by signals, or they have no bus",
                                    device_labels_[device]));
        }
      }

      if (!checks_.passed())
      {
        for (const auto& error : checks_.errors())
        {
          Log::error() << error << "\n";
        }
        throw std::runtime_error("PartitionPlacement::apply failed");
      }

      // External buses stay in the top system
      for (IdxT bus = 0; bus < data_.bus.size(); ++bus)
      {
        if (bus_internal_leaf_[bus] == INVALID_INDEX<IdxT>)
        {
          top.bus.push_back(data_.bus[bus]);
        }
      }

      // Signals of external buses and unused signals stay in the top system
      for (IdxT signal = 0; signal < data_.signal.size(); ++signal)
      {
        const IdxT signal_id   = data_.signal[signal].signal_id;
        IdxT       signal_leaf = INVALID_INDEX<IdxT>;
        if (signal_bus_.contains(signal_id))
        {
          const IdxT bus = signal_bus_.at(signal_id);
          signal_leaf    = bus_internal_leaf_[bus];
        }
        else if (signal_group_.contains(signal_id))
        {
          const IdxT group = signal_group_.at(signal_id);
          signal_leaf      = group_leaf_[group];
        }
        if (signal_leaf == INVALID_INDEX<IdxT>)
        {
          top.signal.push_back(data_.signal[signal]);
        }
      }

      return top;
    }

    /**
     * @brief Build the subsystems of `partition`, or fill `system` if it is a leaf
     *
     * Leaves are numbered in file order, so each signal group goes to the
     * first leaf listing all its buses.
     */
    template <typename real_type, typename index_type>
    void PartitionPlacement<real_type, index_type>::buildSystem(const PartitionDataT& partition,
                                                                DataT&                system)
    {
      if (!partition.partitions.empty())
      {
        // A system with subsystems cannot hold devices on its own buses.
        checks_.check(partition.internal.empty() && partition.external.empty(),
                      std::format("{} has child partitions and bus lists, which is not supported",
                                  partition.name));

        // Subsystem IDs extend the parent ID, as in "0", "1", then "0.0", "0.1"
        std::string prefix;
        if (system.level != 0)
        {
          prefix = system.id + ".";
        }

        system.subsystems.resize(partition.partitions.size());

        // One subsystem per child partition
        for (IdxT child = 0; child < partition.partitions.size(); ++child)
        {
          DataT& subsystem    = system.subsystems[child];
          subsystem.freq_base = system.freq_base;
          subsystem.va_base   = system.va_base;
          subsystem.id        = prefix + std::to_string(child);
          subsystem.level     = system.level + 1;
          buildSystem(partition.partitions[child], subsystem);
        }
      }
      else
      {
        const IdxT leaf = num_leaves_;
        ++num_leaves_;
        addBusList(partition.internal, partition.name, true, leaf);
        addBusList(partition.external, partition.name, false, leaf);

        // Internal buses of this leaf
        for (IdxT bus = 0; bus < data_.bus.size(); ++bus)
        {
          if (bus_internal_leaf_[bus] == leaf)
          {
            system.bus.push_back(data_.bus[bus]);
          }
        }

        // Signal groups not placed yet that have all their buses in this leaf
        for (IdxT device = 0; device < device_group_.size(); ++device)
        {
          if (device_group_[device] == device && group_leaf_[device] == INVALID_INDEX<IdxT>)
          {
            bool has_bus          = false;
            bool all_buses_listed = true;

            // Buses of each device in the group
            for (IdxT member = device;
                 all_buses_listed && member != INVALID_INDEX<IdxT>;
                 member = next_in_group_[member])
            {
              const std::vector<IdxT>& buses = device_buses_[member];
              for (IdxT entry = 0; all_buses_listed && entry < buses.size(); ++entry)
              {
                const IdxT bus   = buses[entry];
                has_bus          = true;
                all_buses_listed = bus_last_leaf_[bus] == leaf;
              }
            }
            if (has_bus && all_buses_listed)
            {
              group_leaf_[device] = leaf;
            }
          }
        }

        CopyLeafDevices<IdxT> copy_devices{device_group_, group_leaf_, leaf, 0};
        DataT::forEachDeviceList(copy_devices, data_, system);

        // Signals of the internal buses and of the groups placed here
        for (IdxT signal = 0; signal < data_.signal.size(); ++signal)
        {
          const IdxT signal_id   = data_.signal[signal].signal_id;
          IdxT       signal_leaf = INVALID_INDEX<IdxT>;
          if (signal_bus_.contains(signal_id))
          {
            const IdxT bus = signal_bus_.at(signal_id);
            signal_leaf    = bus_internal_leaf_[bus];
          }
          else if (signal_group_.contains(signal_id))
          {
            const IdxT group = signal_group_.at(signal_id);
            signal_leaf      = group_leaf_[group];
          }
          if (signal_leaf == leaf)
          {
            system.signal.push_back(data_.signal[signal]);
          }
        }
      }
    }

    /**
     * @brief Record the buses of one list of a leaf partition.
     *
     * `internal` tells whether `buses` is the internal or the external list of
     * leaf `leaf`. A bus internal to a leaf may appear in no other list, and an
     * external bus at most once per leaf. Leaves are listed one after the
     * other, so a bus repeated within `leaf` already has `leaf` as its last
     * leaf.
     */
    template <typename real_type, typename index_type>
    void PartitionPlacement<real_type, index_type>::addBusList(const std::vector<IdxT>& buses,
                                                               const std::string&       name,
                                                               bool                     internal,
                                                               IdxT                     leaf)
    {
      for (const IdxT bus_id : buses)
      {
        checks_.check(bus_index_.contains(bus_id),
                      std::format("{} lists unknown bus {}", name, bus_id));
        if (bus_index_.contains(bus_id))
        {
          const IdxT bus      = bus_index_.at(bus_id);
          const bool repeated = bus_internal_leaf_[bus] != INVALID_INDEX<IdxT>
                                || (internal && bus_last_leaf_[bus] != INVALID_INDEX<IdxT>)
                                || bus_last_leaf_[bus] == leaf;
          checks_.check(!repeated,
                        std::format("{} repeats bus {}, or the bus is internal to two partitions"
                                    " or both internal and external",
                                    name,
                                    bus_id));

          bus_last_leaf_[bus] = leaf;
          if (internal)
          {
            bus_internal_leaf_[bus] = leaf;
          }
        }
      }
    }

    template class PartitionPlacement<double, size_t>;
  } // namespace PhasorDynamics
} // namespace GridKit
