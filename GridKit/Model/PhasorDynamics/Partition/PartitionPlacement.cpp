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

      addAllDevices();

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
        for (IdxT use = 0; use < signal_use_id_.size(); ++use)
        {
          const IdxT device    = signal_use_device_[use];
          const IdxT signal_id = signal_use_id_[use];
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

      group_leaf_.assign(num_devices, INVALID_INDEX<IdxT>);
    }

    /**
     * @brief Add the devices of all types, in the order of fillSystem().
     *
     * Devices are numbered in this order.
     */
    template <typename real_type, typename index_type>
    void PartitionPlacement<real_type, index_type>::addAllDevices()
    {
      addDevices(data_.branch);
      addDevices(data_.regca);
      addDevices(data_.reecb);
      addDevices(data_.repca);
      addDevices(data_.genrou);
      addDevices(data_.gensal);
      addDevices(data_.genclassical);
      addDevices(data_.loadz);
      addDevices(data_.loadzip);
      addDevices(data_.gov);
      addDevices(data_.esdc1a);
      addDevices(data_.gastpti);
      addDevices(data_.hygov);
      addDevices(data_.exciter);
      addDevices(data_.sexspti);
      addDevices(data_.stabilizer);
      addDevices(data_.constant_source);
      addDevices(data_.function_source);
    }

    /**
     * @brief Record the label, buses and signals of each device
     *
     * Buses are stored as case bus indices. Each signal port of a device
     * appends one entry to `signal_use_device_` (the device index) and
     * `signal_use_id_` (the signal ID).
     *
     * @todo Place devices using a bus signal in the leaf where the bus is internal.
     */
    template <typename real_type, typename index_type>
    template <class DeviceDataT>
    void PartitionPlacement<real_type, index_type>::addDevices(
        const std::vector<DeviceDataT>& devices)
    {
      for (const auto& device_data : devices)
      {
        const IdxT device = device_labels_.size();
        device_labels_.push_back(device_data.label());
        device_buses_.emplace_back();

        // Case bus indices of the device's buses
        for (const auto& port : device_data.buses)
        {
          checks_.check(bus_index_.contains(port.second),
                        std::format("device {} uses unknown bus {}",
                                    device_data.label(),
                                    port.second));
          if (bus_index_.contains(port.second))
          {
            device_buses_[device].push_back(bus_index_.at(port.second));
          }
        }

        // A signal of a bus stays in the bus's system, and subsystems do not
        // see their parent's signals. Devices using one are rejected for now.
        for (const auto& port : device_data.signal_inputs)
        {
          signal_use_device_.push_back(device);
          signal_use_id_.push_back(port.second);
          checks_.check(!signal_bus_.contains(port.second),
                        std::format("device {} uses signal {} of a bus",
                                    device_data.label(),
                                    port.second));
        }
        for (const auto& port : device_data.signal_outputs)
        {
          signal_use_device_.push_back(device);
          signal_use_id_.push_back(port.second);
          checks_.check(!signal_bus_.contains(port.second),
                        std::format("device {} uses signal {} of a bus",
                                    device_data.label(),
                                    port.second));
        }
      }
    }

    /**
     * @brief Build the system tree of `partition`, check the placement and fill the systems
     *
     * Every case bus must be listed and every signal group placed in a leaf.
     *
     * @throws std::runtime_error for a badly formed partition
     */
    template <typename real_type, typename index_type>
    SystemModelData<real_type, index_type> PartitionPlacement<real_type, index_type>::place(
        const PartitionDataT& partition)
    {
      // Copy, then clear the buses, signals and devices. Case metadata and
      // monitor sinks stay in the top system.
      DataT top = data_;
      top.bus.clear();
      top.signal.clear();
      top.branch.clear();
      top.regca.clear();
      top.reecb.clear();
      top.repca.clear();
      top.genrou.clear();
      top.gensal.clear();
      top.genclassical.clear();
      top.loadz.clear();
      top.loadzip.clear();
      top.gov.clear();
      top.esdc1a.clear();
      top.gastpti.clear();
      top.hygov.clear();
      top.exciter.clear();
      top.sexspti.clear();
      top.stabilizer.clear();
      top.constant_source.clear();
      top.function_source.clear();

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

      // A signal goes with its bus, or else with the devices using it.
      // Unused signals stay in the top system.
      signal_leaf_.assign(data_.signal.size(), INVALID_INDEX<IdxT>);
      for (IdxT signal = 0; signal < data_.signal.size(); ++signal)
      {
        const IdxT signal_id = data_.signal[signal].signal_id;
        if (signal_bus_.contains(signal_id))
        {
          const IdxT bus       = signal_bus_.at(signal_id);
          signal_leaf_[signal] = bus_internal_leaf_[bus];
        }
        else if (signal_group_.contains(signal_id))
        {
          const IdxT group     = signal_group_.at(signal_id);
          signal_leaf_[signal] = group_leaf_[group];
        }
      }

      // Leaves before the top system, so a root leaf lists its own buses first
      for (IdxT leaf = 0; leaf < leaf_systems_.size(); ++leaf)
      {
        DataT& system = *leaf_systems_[leaf];
        fillSystem(leaf, system);
      }
      fillSystem(INVALID_INDEX<IdxT>, top);

      return top;
    }

    /**
     * @brief Build the subsystems of `partition`, or record its bus lists if it is a leaf
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

        // The subsystem vector is not resized afterwards, so leaf_systems_
        // pointers into it stay valid.
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
        const IdxT leaf = leaf_systems_.size();
        leaf_systems_.push_back(&system);
        addBusList(partition.internal, partition.name, true, leaf);
        addBusList(partition.external, partition.name, false, leaf);

        // Whether each signal group has a bus, and all its buses in this leaf
        std::vector<bool> group_has_bus(device_group_.size(), false);
        std::vector<bool> group_all_listed(device_group_.size(), true);
        for (IdxT device = 0; device < device_group_.size(); ++device)
        {
          const IdxT               group = device_group_[device];
          const std::vector<IdxT>& buses = device_buses_[device];
          for (IdxT entry = 0; entry < buses.size(); ++entry)
          {
            const IdxT bus          = buses[entry];
            group_has_bus[group]    = true;
            group_all_listed[group] = group_all_listed[group] && bus_last_leaf_[bus] == leaf;
          }
        }

        // Signal groups not placed yet that have all their buses in this leaf
        for (IdxT group = 0; group < device_group_.size(); ++group)
        {
          if (device_group_[group] == group && group_leaf_[group] == INVALID_INDEX<IdxT>
              && group_has_bus[group] && group_all_listed[group])
          {
            group_leaf_[group] = leaf;
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

    /**
     * @brief Copy the buses, signals and devices of `leaf` from the flat case into `system`
     *
     * `leaf` is INVALID_INDEX for the top system, which gets the external
     * buses and the signals of no leaf. Device types are in the order of
     * addAllDevices().
     */
    template <typename real_type, typename index_type>
    void PartitionPlacement<real_type, index_type>::fillSystem(IdxT leaf, DataT& system) const
    {
      for (IdxT bus = 0; bus < data_.bus.size(); ++bus)
      {
        if (bus_internal_leaf_[bus] == leaf)
        {
          system.bus.push_back(data_.bus[bus]);
        }
      }

      for (IdxT signal = 0; signal < data_.signal.size(); ++signal)
      {
        if (signal_leaf_[signal] == leaf)
        {
          system.signal.push_back(data_.signal[signal]);
        }
      }

      IdxT device = 0;
      copyDevices(data_.branch, leaf, system.branch, device);
      copyDevices(data_.regca, leaf, system.regca, device);
      copyDevices(data_.reecb, leaf, system.reecb, device);
      copyDevices(data_.repca, leaf, system.repca, device);
      copyDevices(data_.genrou, leaf, system.genrou, device);
      copyDevices(data_.gensal, leaf, system.gensal, device);
      copyDevices(data_.genclassical, leaf, system.genclassical, device);
      copyDevices(data_.loadz, leaf, system.loadz, device);
      copyDevices(data_.loadzip, leaf, system.loadzip, device);
      copyDevices(data_.gov, leaf, system.gov, device);
      copyDevices(data_.esdc1a, leaf, system.esdc1a, device);
      copyDevices(data_.gastpti, leaf, system.gastpti, device);
      copyDevices(data_.hygov, leaf, system.hygov, device);
      copyDevices(data_.exciter, leaf, system.exciter, device);
      copyDevices(data_.sexspti, leaf, system.sexspti, device);
      copyDevices(data_.stabilizer, leaf, system.stabilizer, device);
      copyDevices(data_.constant_source, leaf, system.constant_source, device);
      copyDevices(data_.function_source, leaf, system.function_source, device);

      // A type missing above would leave devices uncopied.
      if (device != device_labels_.size())
      {
        Log::error() << "Device types of fillSystem and addAllDevices differ\n";
        throw std::runtime_error("PartitionPlacement::apply failed");
      }
    }

    /**
     * @brief Copy the devices of one type whose signal group is placed in `leaf`
     *
     * `device` numbers the devices of all types, as in addAllDevices().
     */
    template <typename real_type, typename index_type>
    template <class DeviceDataT>
    void PartitionPlacement<real_type, index_type>::copyDevices(
        const std::vector<DeviceDataT>& devices,
        IdxT                            leaf,
        std::vector<DeviceDataT>&       system_devices,
        IdxT&                           device) const
    {
      for (const auto& device_data : devices)
      {
        // Device types listed in another order than in addAllDevices() would
        // give devices the groups of others.
        if (device >= device_labels_.size() || device_labels_[device] != device_data.label())
        {
          Log::error() << "Device " << device_data.label()
                       << " is out of order in fillSystem compared to addAllDevices\n";
          throw std::runtime_error("PartitionPlacement::apply failed");
        }

        const IdxT group = device_group_[device];
        if (group_leaf_[group] == leaf)
        {
          system_devices.push_back(device_data);
        }
        ++device;
      }
    }

    template class PartitionPlacement<double, size_t>;
  } // namespace PhasorDynamics
} // namespace GridKit
