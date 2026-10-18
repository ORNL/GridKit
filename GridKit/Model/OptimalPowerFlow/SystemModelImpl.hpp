/**
 * @file SystemModelImpl.hpp
 * @brief Definition of the optimal power flow system model.
 */

#pragma once

#include <format>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_set>
#include <utility>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/OptimalPowerFlow/Branch/Branch.hpp>
#include <GridKit/Model/OptimalPowerFlow/Bus/Bus.hpp>
#include <GridKit/Model/OptimalPowerFlow/ComponentImpl.hpp>
#include <GridKit/Model/OptimalPowerFlow/Generator/Generator.hpp>
#include <GridKit/Model/OptimalPowerFlow/Load/Load.hpp>
#include <GridKit/Model/OptimalPowerFlow/SystemModel.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /**
     * @param[in] data - Network, demand, limits, and costs
     *
     * @throw std::invalid_argument for a duplicate bus, a parameter that
     * cannot be read, a port without a bus, or a bus that does not exist
     */
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::SystemModel(const SystemModelData<RealT, IdxT>& data)
    {
      for (const auto& bus : data.bus)
      {
        addBus(bus);
      }

      for (const auto& branch : data.branch)
      {
        branches_.push_back(&addDevice<BranchT>(branch));
      }

      for (const auto& generator : data.generator)
      {
        addDevice<Generator<ScalarT, IdxT>>(generator);
      }

      for (const auto& load : data.load)
      {
        addDevice<Load<ScalarT, IdxT>>(load);
      }
    }

    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::SystemModel(std::vector<std::unique_ptr<ComponentT>> components)
      : components_(std::move(components))
    {
      for (const auto& component : components_)
      {
        if (!component)
        {
          throw std::invalid_argument("OptimalPowerFlow::SystemModel: null component");
        }
      }
    }

    /**
     * @brief Every problem with the data, including an island without a
     * reference bus
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks SystemModel<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks;

      std::unordered_set<const ComponentT*> components;
      for (const auto& component : components_)
      {
        components.insert(component.get());
      }

      for (const auto& [number, bus] : buses_)
      {
        checks.check(number != INVALID_INDEX<IdxT>,
                     std::format("{}: invalid bus number", bus->id()));
      }

      for (const IdxT number : unreferencedIslands())
      {
        checks.fail(std::format("OptimalPowerFlow::SystemModel: the island of bus {} has no reference bus", number));
      }

      std::set<std::string> ids;
      for (const auto& component : components_)
      {
        checks.check(ids.insert(component->id()).second, std::format("{}: duplicate id", component->id()));

        const auto component_checks = component->verify();
        for (const auto& error : component_checks.errors())
        {
          checks.fail(std::format("{}: {}", component->id(), error));
        }

        const auto& ports = component->ports();
        for (IdxT port = 0; port < ports.size(); ++port)
        {
          checks.check(components.contains(ports[port].component),
                       std::format("{}: port {} must connect to a component in this system", component->id(), port));
        }
      }

      return checks;
    }

    /**
     * @brief Verify the components, bind them to the system vectors, and
     * build the derivative patterns
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::allocate()
    {
      const auto checks = verify();
      for (const auto& error : checks.errors())
      {
        Log::error() << error << '\n';
      }
      if (!checks.passed())
      {
        return 1;
      }

      if (compose() != 0)
      {
        return 1;
      }

      return allocateDerivatives();
    }

    /**
     * @brief Bounds and the flat start
     *
     * Every variable starts at zero and every bound unbounded. Each component
     * then sets its bounds and any nonzero initial value.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::initialize()
    {
      if (x_.setToZero() != 0
          || x_lower_.setToConst(-UNBOUNDED<RealT>) != 0
          || x_upper_.setToConst(UNBOUNDED<RealT>) != 0
          || d_lower_.setToConst(-UNBOUNDED<RealT>) != 0
          || d_upper_.setToConst(UNBOUNDED<RealT>) != 0)
      {
        return 1;
      }

      for (auto& component : components_)
      {
        const int status = component->initialize();
        if (status != 0)
        {
          return status;
        }
      }

      // For DependencyTracking::Variable, set variable numbers
      if constexpr (std::is_same_v<ScalarT, DependencyTracking::Variable>)
      {
        ScalarT* x = x_.getData();
        for (IdxT j = 0; j < x_.getSize(); ++j)
        {
          x[j].setVariableNumber(j);
        }
      }

      x_.setDataUpdated();

      return 0;
    }

    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateObjective()
    {
      f_ = ZERO<RealT>;
      for (auto& component : components_)
      {
        const int status = component->evaluateObjective();
        if (status != 0)
        {
          return status;
        }
        f_ += component->f();
      }

      return 0;
    }

    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateEqualities()
    {
      c_.setToZero();
      for (auto& component : components_)
      {
        const int status = component->evaluateEqualities();
        if (status != 0)
        {
          return status;
        }
      }
      c_.setDataUpdated();

      return 0;
    }

    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateInequalities()
    {
      d_.setToZero();
      for (auto& component : components_)
      {
        const int status = component->evaluateInequalities();
        if (status != 0)
        {
          return status;
        }
      }
      d_.setDataUpdated();

      return 0;
    }

    /**
     * @brief `ModelT` from `data`, with a parameter error naming `id`
     */
    template <typename scalar_type, typename index_type>
    template <typename ModelT, typename DataT>
    std::unique_ptr<ModelT> SystemModel<scalar_type, index_type>::build(const DataT& data, std::string_view id)
    {
      try
      {
        return std::make_unique<ModelT>(data);
      }
      catch (const std::invalid_argument& error)
      {
        throw std::invalid_argument(std::format("{}: {}", id, error.what()));
      }
    }

    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::addBus(const BusData<RealT, IdxT>& data)
    {
      auto bus = build<BusT>(data, data.id());
      if (!buses_.emplace(data.number, bus.get()).second)
      {
        throw std::invalid_argument(std::format("OptimalPowerFlow::SystemModel: duplicate bus {}", data.number));
      }
      components_.push_back(std::move(bus));
    }

    /**
     * @brief Device from `data`, connected at each port to the bus it names
     */
    template <typename scalar_type, typename index_type>
    template <typename ModelT, typename DataT>
    ModelT& SystemModel<scalar_type, index_type>::addDevice(const DataT& data)
    {
      using Ports = typename DataT::Ports;
      static_assert(ModelT::NUM_PORT_VARIABLES.size() == Utilities::enum_size<Ports>());

      auto device = build<ModelT>(data, data.id);
      for (const Ports port : Utilities::enum_values<Ports>())
      {
        const auto number = data.ports.find(port);
        if (number == data.ports.end())
        {
          throw std::invalid_argument(std::format("{}: port '{}' has no bus", data.id, Utilities::enum_name(port)));
        }

        const auto bus = buses_.find(number->second);
        if (bus == buses_.end())
        {
          throw std::invalid_argument(
              std::format("{}: port '{}' names unknown bus {}", data.id, Utilities::enum_name(port), number->second));
        }

        if (device->connect(static_cast<IdxT>(port), *bus->second) != 0)
        {
          throw std::invalid_argument(
              std::format("{}: port '{}' does not fit bus {}", data.id, Utilities::enum_name(port), number->second));
        }
      }

      ModelT& added = *device;
      components_.push_back(std::move(device));
      return added;
    }

    /**
     * @brief Size the vectors for the components and bind each to the next
     * blocks
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::compose()
    {
      Block size;
      for (const auto& component : components_)
      {
        size += component->size_;
      }

      if (x_.resize(size.variables) != 0
          || x_lower_.resize(size.variables) != 0
          || x_upper_.resize(size.variables) != 0
          || gradient_.resize(size.variables) != 0
          || c_.resize(size.equalities) != 0
          || d_.resize(size.inequalities) != 0
          || d_lower_.resize(size.inequalities) != 0
          || d_upper_.resize(size.inequalities) != 0)
      {
        return 1;
      }

      Block offset;
      for (const auto& component : components_)
      {
        if (bind(*component, offset) != 0)
        {
          return 1;
        }
        offset += component->size_;
      }
      return 0;
    }

    /// View the internal blocks of the system vectors from `offset`
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::bind(ComponentT& component, const Block& offset)
    {
      component.offset_ = offset;
      const Block& size = component.size_;

      const auto view = [](auto& local, auto& global, IdxT start, IdxT count)
      { return local.setData(global.getData() + start, count) != 0; };

      if (view(component.x_, x_, offset.variables, size.variables)
          || view(component.x_lower_, x_lower_, offset.variables, size.variables)
          || view(component.x_upper_, x_upper_, offset.variables, size.variables)
          || view(component.gradient_, gradient_, offset.variables, size.variables)
          || view(component.c_, c_, offset.equalities, size.equalities)
          || view(component.d_, d_, offset.inequalities, size.inequalities)
          || view(component.d_lower_, d_lower_, offset.inequalities, size.inequalities)
          || view(component.d_upper_, d_upper_, offset.inequalities, size.inequalities))
      {
        return 1;
      }
      component.bound_ = true;
      return 0;
    }

    /**
     * @brief Lowest bus number of each island without a reference bus
     *
     * In-service branches connect the buses of their ports.
     */
    template <typename scalar_type, typename index_type>
    std::vector<index_type> SystemModel<scalar_type, index_type>::unreferencedIslands() const
    {
      // Each bus points toward the root bus of its island
      std::map<const ComponentT*, const ComponentT*> parent;
      for (const auto& [number, bus] : buses_)
      {
        parent[bus] = bus;
      }

      const auto root = [&parent](const ComponentT* bus)
      {
        while (parent.at(bus) != bus)
        {
          bus = parent.at(bus);
        }
        return bus;
      };

      for (const auto* branch : branches_)
      {
        const auto& ports = branch->ports();
        if (!branch->inService() || ports[0].component == nullptr || ports[1].component == nullptr)
        {
          continue;
        }
        parent[root(ports[1].component)] = root(ports[0].component);
      }

      std::set<const ComponentT*> referenced;
      for (const auto& [number, bus] : buses_)
      {
        if (bus->reference())
        {
          referenced.insert(root(bus));
        }
      }

      std::vector<IdxT>           unreferenced;
      std::set<const ComponentT*> reported;
      for (const auto& [number, bus] : buses_)
      {
        const auto* island = root(bus);
        if (!referenced.contains(island) && reported.insert(island).second)
        {
          unreferenced.push_back(number);
        }
      }
      return unreferenced;
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
