/**
 * @file SystemModelImpl.hpp
 * @brief Definition of the optimal power flow system model.
 */

#pragma once

#include <algorithm>
#include <format>
#include <map>
#include <set>
#include <stdexcept>
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

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /**
     * @param[in] data - Network, demand, limits, and costs
     */
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::SystemModel(const SystemModelData<RealT, IdxT>& data)
    {
      for (const auto& bus_data : data.bus)
      {
        auto bus = std::make_unique<BusT>(bus_data);
        if (!buses_.emplace(bus_data.number, bus.get()).second)
        {
          throw std::invalid_argument(std::format("OptimalPowerFlow::SystemModel: duplicate bus {}", bus_data.number));
        }
        components_.push_back(std::move(bus));
      }

      const auto add = [this](auto component, const auto& component_data)
      {
        const auto numbers = portNumbers(component_data);
        for (IdxT port = 0; port < numbers.size(); ++port)
        {
          const auto bus = buses_.find(numbers[port]);
          if (numbers[port] != INVALID_INDEX<IdxT> && bus != buses_.end())
          {
            const int status = component->connect(port, *bus->second);
            if (status != 0)
            {
              throw std::invalid_argument(std::format("{}: incompatible port {}", component->id(), port));
            }
          }
        }
        components_.push_back(std::move(component));
      };

      for (const auto& branch_data : data.branch)
      {
        auto branch = std::make_unique<BranchT>(branch_data);
        branches_.push_back(branch.get());
        add(std::move(branch), branch_data);
      }

      for (const auto& generator : data.generator)
      {
        add(std::make_unique<Generator<ScalarT, IdxT>>(generator), generator);
      }

      for (const auto& load : data.load)
      {
        add(std::make_unique<Load<ScalarT, IdxT>>(load), load);
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

    /// View the internal blocks of the system vectors from `offset`
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::bind(ComponentT& component, const Block& offset)
    {
      component.offset_ = offset;
      const Block& size = component.size_;

      const auto view = [](auto& local, auto& global, IdxT start, IdxT size)
      { return local.setData(global.getData() + start, size) != 0; };

      if (view(component.x_, this->x_, offset.variables, size.variables)
          || view(component.x_lower_, this->x_lower_, offset.variables, size.variables)
          || view(component.x_upper_, this->x_upper_, offset.variables, size.variables)
          || view(component.gradient_, this->gradient_, offset.variables, size.variables)
          || view(component.c_, this->c_, offset.equalities, size.equalities)
          || view(component.d_, this->d_, offset.inequalities, size.inequalities)
          || view(component.d_lower_, this->d_lower_, offset.inequalities, size.inequalities)
          || view(component.d_upper_, this->d_upper_, offset.inequalities, size.inequalities))
      {
        return 1;
      }
      component.bound_ = true;
      return 0;
    }

    /**
     * @brief Size the vectors for `components`, bind each to the next blocks,
     * and allocate it
     *
     * Limits start unlimited, and each component then sets its limits.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::compose()
    {
      Block size;
      for (const auto& component : components_)
      {
        size += component->size_;
      }

      if (this->x_.resize(size.variables) != 0
          || this->x_lower_.resize(size.variables) != 0
          || this->x_upper_.resize(size.variables) != 0
          || this->gradient_.resize(size.variables) != 0
          || this->c_.resize(size.equalities) != 0
          || this->d_.resize(size.inequalities) != 0
          || this->d_lower_.resize(size.inequalities) != 0
          || this->d_upper_.resize(size.inequalities) != 0)
      {
        return 1;
      }

      if (this->x_lower_.setToConst(-UNBOUNDED<RealT>) != 0
          || this->x_upper_.setToConst(UNBOUNDED<RealT>) != 0
          || this->d_lower_.setToConst(-UNBOUNDED<RealT>) != 0
          || this->d_upper_.setToConst(UNBOUNDED<RealT>) != 0)
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
      for (const auto& component : components_)
      {
        const int status = component->allocate();
        if (status != 0)
        {
          return status;
        }
      }
      return 0;
    }

    /**
     * @brief Compose the components and set the
     * starting point and the derivative patterns
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

      const int status = initialize();
      if (status != 0)
      {
        return status;
      }

      return allocateDerivatives();
    }

    /**
     * @brief Starting point
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::initialize()
    {
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

      for (const auto& component : components_)
      {
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
