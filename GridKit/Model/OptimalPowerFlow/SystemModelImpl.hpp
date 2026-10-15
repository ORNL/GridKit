/**
 * @file SystemModelImpl.hpp
 * @brief Definition of the optimal power flow system model.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <format>
#include <map>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/OptimalPowerFlow/Branch/Branch.hpp>
#include <GridKit/Model/OptimalPowerFlow/Bus/Bus.hpp>
#include <GridKit/Model/OptimalPowerFlow/Generator/Generator.hpp>
#include <GridKit/Model/OptimalPowerFlow/Load/Load.hpp>
#include <GridKit/Model/OptimalPowerFlow/Shunt/Shunt.hpp>
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

      for (const auto& branch : data.branch)
      {
        components_.push_back(std::make_unique<Branch<ScalarT, IdxT>>(branch));
      }

      for (const auto& generator : data.generator)
      {
        components_.push_back(std::make_unique<Generator<ScalarT, IdxT>>(generator));
      }

      for (const auto& load : data.load)
      {
        components_.push_back(std::make_unique<Load<ScalarT, IdxT>>(load));
      }

      for (const auto& shunt : data.shunt)
      {
        components_.push_back(std::make_unique<Shunt<ScalarT, IdxT>>(shunt));
      }
    }

    /**
     * @brief Number the variables and the limited constraint functions,
     * connect the terminals, copy the limits, and set the starting point and
     * the derivative patterns
     *
     * A constraint function without a finite limit gets no row.
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

      IdxT n = 0;
      IdxT m = 0;
      for (auto& component : components_)
      {
        auto& variables = component->variableIndices();
        for (IdxT j = 0; j < component->sizeInternalVariables(); ++j)
        {
          variables[j] = n++;
        }

        auto& rows = component->constraintIndices();
        for (IdxT i = 0; i < component->sizeInternalConstraints(); ++i)
        {
          rows[i] = INVALID_INDEX<IdxT>;
          if (std::isfinite(component->constraintLower()[i]) || std::isfinite(component->constraintUpper()[i]))
          {
            rows[i] = m++;
          }
        }
      }

      for (auto& component : components_)
      {
        if (connect(*component) != 0)
        {
          return 1;
        }
      }

      x_.resize(n);
      x_lower_.resize(n);
      x_upper_.resize(n);
      gradient_.resize(n);
      g_.resize(m);
      g_lower_.resize(m);
      g_upper_.resize(m);

      RealT* x_lower = x_lower_.getData();
      RealT* x_upper = x_upper_.getData();
      RealT* g_lower = g_lower_.getData();
      RealT* g_upper = g_upper_.getData();
      for (const auto& component : components_)
      {
        const auto& variables = component->variableIndices();
        for (IdxT j = 0; j < component->sizeInternalVariables(); ++j)
        {
          x_lower[variables[j]] = component->variableLower()[j];
          x_upper[variables[j]] = component->variableUpper()[j];
        }

        const auto& rows = component->constraintIndices();
        for (IdxT i = 0; i < component->sizeInternalConstraints(); ++i)
        {
          if (rows[i] != INVALID_INDEX<IdxT>)
          {
            g_lower[rows[i]] = component->constraintLower()[i];
            g_upper[rows[i]] = component->constraintUpper()[i];
          }
        }
      }
      x_lower_.setDataUpdated();
      x_upper_.setDataUpdated();
      g_lower_.setDataUpdated();
      g_upper_.setDataUpdated();

      if (initialize() != 0)
      {
        return 1;
      }

      return allocateDerivatives();
    }

    /**
     * @brief Starting point
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::initialize()
    {
      ScalarT* x = x_.getData();

      int ret = 0;
      for (auto& component : components_)
      {
        ret += component->initialize(x);
      }

      // For DependencyTracking::Variable, set variable numbers
      if constexpr (std::is_same_v<ScalarT, DependencyTracking::Variable>)
      {
        for (IdxT j = 0; j < x_.getSize(); ++j)
        {
          x[j].setVariableNumber(j);
        }
      }

      x_.setDataUpdated();

      return ret;
    }

    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateObjective()
    {
      const ScalarT* x = x_.getData();

      f_ = ZERO<RealT>;
      for (auto& component : components_)
      {
        component->evaluateObjective(x, f_);
      }

      return 0;
    }

    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateConstraints()
    {
      const ScalarT* x = x_.getData();

      g_.setToZero();
      ScalarT* g = g_.getData();
      for (auto& component : components_)
      {
        component->evaluateConstraints(x, g);
      }
      g_.setDataUpdated();

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

      for (const IdxT number : unreferencedIslands())
      {
        checks.fail(std::format("OptimalPowerFlow::SystemModel: the island of bus {} has no reference bus", number));
      }

      for (const auto& component : components_)
      {
        for (const auto& error : component->verify().errors())
        {
          checks.fail(std::format("{}: {}", component->id(), error));
        }
      }

      return checks;
    }

    /**
     * @brief Lowest bus number of each island without a reference bus
     *
     * In-service devices connect the buses of their terminals.
     */
    template <typename scalar_type, typename index_type>
    std::vector<index_type> SystemModel<scalar_type, index_type>::unreferencedIslands() const
    {
      // Each bus points toward the root bus of its island
      std::map<IdxT, IdxT> parent;
      for (const auto& [number, bus] : buses_)
      {
        parent[number] = number;
      }

      const auto root = [&parent](IdxT number)
      {
        while (parent.at(number) != number)
        {
          number = parent.at(number);
        }
        return number;
      };

      for (const auto& component : components_)
      {
        const auto& terminals = component->terminals();
        const bool  connected = std::ranges::all_of(terminals, [this](IdxT number)
                                                   { return buses_.contains(number); });
        if (!component->inService() || !connected)
        {
          continue;
        }
        for (IdxT t = 1; t < terminals.size(); ++t)
        {
          parent[root(terminals[t])] = root(terminals[0]);
        }
      }

      std::set<IdxT> referenced;
      for (const auto& [number, bus] : buses_)
      {
        if (bus->reference())
        {
          referenced.insert(root(number));
        }
      }

      std::vector<IdxT> unreferenced;
      std::set<IdxT>    reported;
      for (const auto& [number, bus] : buses_)
      {
        const IdxT island = root(number);
        if (!referenced.contains(island) && reported.insert(island).second)
        {
          unreferenced.push_back(number);
        }
      }
      return unreferenced;
    }

    /**
     * @brief Terminal slots take the voltage and balance rows of their bus
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::connect(ComponentT& component)
    {
      const auto& terminals = component.terminals();
      for (IdxT t = 0; t < terminals.size(); ++t)
      {
        const auto bus = buses_.find(terminals[t]);
        if (bus == buses_.end())
        {
          Log::error() << component.id() << ": terminal " << t << " has no bus\n";
          return 1;
        }

        const IdxT variable   = component.sizeInternalVariables() + ComponentT::TERMINAL_SIZE * t;
        const IdxT constraint = component.sizeInternalConstraints() + ComponentT::TERMINAL_SIZE * t;
        for (IdxT k = 0; k < ComponentT::TERMINAL_SIZE; ++k)
        {
          component.variableIndices()[variable + k]     = bus->second->variableIndices()[k];
          component.constraintIndices()[constraint + k] = bus->second->constraintIndices()[k];
        }
      }
      return 0;
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
