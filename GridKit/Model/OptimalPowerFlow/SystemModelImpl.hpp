/**
 * @file SystemModelImpl.hpp
 * @brief Definition of the optimal power flow system model.
 */

#pragma once

#include <format>
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
        addDevice<Branch<ScalarT, IdxT>>(branch);
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
     * @brief Every problem with the data
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
     * @brief Verify the components, bind them to the system vectors, tag the
     * inequalities, and build the derivative patterns
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

      if (compose() != 0 || tagInequality() != 0)
      {
        return 1;
      }

      return allocateDerivatives();
    }

    /**
     * @brief Each component tags its own constraints
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::tagInequality()
    {
      for (auto& component : components_)
      {
        const int status = component->tagInequality();
        if (status != 0)
        {
          return status;
        }

        const auto& tag = component->tag();
        for (IdxT i = 0; i < tag.size(); ++i)
        {
          tag_[component->offset_.constraints + i] = tag[i];
        }
      }

      return 0;
    }

    /**
     * @brief Bounds and the flat start
     *
     * Every variable and multiplier starts at zero and every variable
     * unbounded, every equality at zero, and every inequality unbounded. Each
     * component then sets its limits and any nonzero initial value.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::initialize()
    {
      if (x_.setToZero() != 0
          || x_lower_.setToConst(-UNBOUNDED<RealT>) != 0
          || x_upper_.setToConst(UNBOUNDED<RealT>) != 0
          || g_lower_.setToZero() != 0
          || g_upper_.setToZero() != 0
          || lambda_.setToZero() != 0
          || z_lower_.setToZero() != 0
          || z_upper_.setToZero() != 0)
      {
        return 1;
      }

      RealT* g_lower = g_lower_.getData();
      RealT* g_upper = g_upper_.getData();
      for (IdxT i = 0; i < g_.getSize(); ++i)
      {
        if (tag_[i])
        {
          g_lower[i] = -UNBOUNDED<RealT>;
          g_upper[i] = UNBOUNDED<RealT>;
        }
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
    int SystemModel<scalar_type, index_type>::evaluateConstraints()
    {
      g_.setToZero();
      for (auto& component : components_)
      {
        const int status = component->evaluateConstraints();
        if (status != 0)
        {
          return status;
        }
      }
      g_.setDataUpdated();

      return 0;
    }

    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::addBus(const BusData<RealT, IdxT>& data)
    {
      auto bus = std::make_unique<BusT>(data);
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
    void SystemModel<scalar_type, index_type>::addDevice(const DataT& data)
    {
      using Ports = typename DataT::Ports;
      static_assert(ModelT::NUM_PORT_VARIABLES.size() == Utilities::enum_size<Ports>());

      auto device = std::make_unique<ModelT>(data);
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

      components_.push_back(std::move(device));
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
          || z_lower_.resize(size.variables) != 0
          || z_upper_.resize(size.variables) != 0
          || g_.resize(size.constraints) != 0
          || g_lower_.resize(size.constraints) != 0
          || g_upper_.resize(size.constraints) != 0
          || lambda_.resize(size.constraints) != 0)
      {
        return 1;
      }
      tag_.assign(size.constraints, false);

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
          || view(component.g_, g_, offset.constraints, size.constraints)
          || view(component.g_lower_, g_lower_, offset.constraints, size.constraints)
          || view(component.g_upper_, g_upper_, offset.constraints, size.constraints)
          || view(component.lambda_, lambda_, offset.constraints, size.constraints))
      {
        return 1;
      }
      component.bound_ = true;
      return 0;
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
