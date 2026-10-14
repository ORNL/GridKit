#pragma once

#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace PartitionDetail
    {
      using DataT = SystemModelData<>;
      using IdxT  = DataT::IdxT;
      using json  = nlohmann::json;

      // Keep the list of component containers in one place for all placement passes.
      inline void forEachComponentList(const auto& visit)
      {
        visit(&DataT::branch);
        visit(&DataT::regca);
        visit(&DataT::reecb);
        visit(&DataT::repca);
        visit(&DataT::genrou);
        visit(&DataT::gensal);
        visit(&DataT::genclassical);
        visit(&DataT::loadz);
        visit(&DataT::loadzip);
        visit(&DataT::gov);
        visit(&DataT::gastpti);
        visit(&DataT::hygov);
        visit(&DataT::exciter);
        visit(&DataT::esdc1a);
        visit(&DataT::sexspti);
        visit(&DataT::stabilizer);
        visit(&DataT::constant_source);
        visit(&DataT::function_source);
      }

      struct Leaf
      {
        DataT*         data;
        std::set<IdxT> buses;
      };

      inline void buildPartition(const json&             node,
                                 DataT&                  data,
                                 DataT&                  root,
                                 std::vector<Leaf>&      leaves,
                                 std::map<IdxT, DataT*>& bus_owner)
      {
        const auto name  = node.at("name").get<std::string>();
        const bool group = node.contains("partitions");
        if (name.empty())
        {
          throw std::runtime_error("Invalid partition object: " + name);
        }
        if (group)
        {
          if (node.contains("internal") || node.contains("external"))
          {
            throw std::runtime_error("Partition group cannot contain bus lists: " + name);
          }
          const auto& children = node.at("partitions");
          if (!children.is_array() || children.empty())
          {
            throw std::runtime_error("Partition group must contain partitions: " + name);
          }
          data.subsystems.resize(children.size());
          for (std::size_t i = 0; i < children.size(); ++i)
          {
            auto& child     = data.subsystems[i];
            child.id        = data.id + "/" + std::to_string(i);
            child.level     = data.level + 1;
            child.freq_base = data.freq_base;
            child.va_base   = data.va_base;
            buildPartition(children[i], child, root, leaves, bus_owner);
          }
          return;
        }

        auto& leaf     = leaves.emplace_back(Leaf{&data, {}});
        auto  addBuses = [&](const json& buses, DataT* owner)
        {
          if (!buses.is_array())
          {
            throw std::runtime_error("Partition bus list must be an array: " + name);
          }
          for (const auto& value : buses)
          {
            if (!value.is_number_integer() || value < 0)
            {
              throw std::runtime_error("Partition bus ID must be a nonnegative integer: " + name);
            }
            const auto id = value.get<IdxT>();
            auto       it = bus_owner.find(id);
            if (it == bus_owner.end() || !leaf.buses.insert(id).second)
            {
              throw std::runtime_error("Unknown or repeated bus " + std::to_string(id) + " in partition " + name);
            }
            if (it->second != nullptr && (owner != &root || it->second != &root))
            {
              throw std::runtime_error("Conflicting ownership of bus " + std::to_string(id));
            }
            it->second = owner;
          }
        };
        addBuses(node.at("internal"), &data);
        addBuses(node.at("external"), &root);
      }
    } // namespace PartitionDetail

    /**
     * @brief Compose a flat case into subsystems using a partition specification.
     *
     * External buses belong to the root. Leaves own their internal buses and
     * devices; signals must stay within one system. Placement is checked before
     * moving data, and case metadata and monitor sinks remain at the root.
     */
    inline void applyPartition(SystemModelData<>& data, const nlohmann::json& partition)
    {
      using DataT = SystemModelData<>;
      using IdxT  = DataT::IdxT;
      using PartitionDetail::forEachComponentList;

      if (!data.subsystems.empty() || !partition.contains("partitions"))
      {
        throw std::runtime_error("Partitioning requires a flat case and a root partition group");
      }

      DataT root;
      root.freq_base = data.freq_base;
      root.va_base   = data.va_base;

      std::vector<PartitionDetail::Leaf> leaves;
      std::map<IdxT, DataT*>             bus_owner;
      for (const auto& bus : data.bus)
      {
        if (!bus_owner.emplace(bus.bus_id, nullptr).second)
        {
          throw std::runtime_error("Duplicate case bus " + std::to_string(bus.bus_id));
        }
      }
      PartitionDetail::buildPartition(partition, root, root, leaves, bus_owner);
      for (const auto& [id, owner] : bus_owner)
      {
        if (owner == nullptr)
        {
          throw std::runtime_error("Bus missing from partition: " + std::to_string(id));
        }
      }

      auto forEachComponent = [&](const auto& visit)
      {
        std::size_t index     = 0;
        auto        visitList = [&](auto member)
        {
          for (auto& component : data.*member)
          {
            visit(component, member, index++);
          }
        };
        forEachComponentList(visitList);
      };

      std::vector<DataT*> component_owner;
      auto                assignComponent = [&](const auto& component, auto, std::size_t)
      {
        DataT* owner = nullptr;
        if (!component.buses.empty())
        {
          for (const auto& leaf : leaves)
          {
            if (std::ranges::all_of(component.buses, [&](const auto& port)
                                    { return leaf.buses.contains(port.second); }))
            {
              owner = leaf.data;
              break;
            }
          }
          if (owner == nullptr)
          {
            throw std::runtime_error("No partition contains all buses of component " + component.disambiguation_string);
          }
        }
        component_owner.push_back(owner);
      };
      forEachComponent(assignComponent);

      std::map<IdxT, DataT*> signal_owner;
      for (const auto& signal : data.signal)
      {
        if (!signal_owner.emplace(signal.signal_id, nullptr).second)
        {
          throw std::runtime_error("Duplicate case signal " + std::to_string(signal.signal_id));
        }
      }
      auto propagateSignalOwnership = [&](const auto& component, DataT*& owner)
      {
        bool changed = false;
        auto ports   = [&](const auto& mappings)
        {
          for (const auto& [port, id] : mappings)
          {
            auto it = signal_owner.find(id);
            if (it == signal_owner.end())
            {
              throw std::runtime_error("Unknown signal " + std::to_string(id));
            }
            auto*& signal = it->second;
            if (owner && signal && owner != signal)
            {
              throw std::runtime_error("Signal crosses partitions: " + std::to_string(id));
            }
            if (!owner && signal)
            {
              owner   = signal;
              changed = true;
            }
            if (owner && !signal)
            {
              signal  = owner;
              changed = true;
            }
          }
        };
        ports(component.signal_inputs);
        ports(component.signal_outputs);
        return changed;
      };
      for (const auto& bus : data.bus)
      {
        auto* owner = bus_owner.at(bus.bus_id);
        propagateSignalOwnership(bus, owner);
      }
      bool changed;
      do
      {
        changed = false;
        forEachComponent([&](const auto& component, auto, std::size_t index)
                         { changed |= propagateSignalOwnership(component, component_owner[index]); });
      } while (changed);

      auto checkComponent = [&](const auto& component, auto, std::size_t index)
      {
        auto* owner = component_owner[index];
        if (owner == nullptr || owner == &root)
        {
          throw std::runtime_error("Component has no leaf partition: " + component.disambiguation_string);
        }
      };
      forEachComponent(checkComponent);

      for (auto& bus : data.bus)
      {
        bus_owner.at(bus.bus_id)->bus.push_back(std::move(bus));
      }
      for (auto& signal : data.signal)
      {
        auto* owner = signal_owner.at(signal.signal_id);
        (owner ? owner : &root)->signal.push_back(std::move(signal));
      }
      forEachComponent([&](auto& component, auto member, std::size_t index)
                       { (component_owner[index]->*member).push_back(std::move(component)); });
      forEachComponentList([&](auto member)
                           { (data.*member).clear(); });
      data.bus        = std::move(root.bus);
      data.signal     = std::move(root.signal);
      data.subsystems = std::move(root.subsystems);
    }
  } // namespace PhasorDynamics
} // namespace GridKit
