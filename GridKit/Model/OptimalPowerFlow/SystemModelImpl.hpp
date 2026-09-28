/**
 * @file SystemModelImpl.hpp
 * @brief Definition of the optimal power flow system model.
 */

#pragma once

#include <cmath>
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
     * @param[in] data - Network, limits, and costs
     * @param[in] state - Starting point, fixed demand, and device settings
     */
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::SystemModel(const SystemModelData<RealT, IdxT>& data,
                                                      const Model::StateData&             state)
      : state_(state)
    {
      for (const auto& bus_data : data.bus)
      {
        auto bus = std::make_unique<BusT>(bus_data);
        if (!buses_.emplace(bus_data.number, bus.get()).second)
        {
          Log::error() << "OptimalPowerFlow::SystemModel: duplicate bus " << bus_data.number << "\n";
          ++data_errors_;
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
     * @brief Number the variables and the bounded constraints, connect the
     * terminals, and set the starting point and the derivative patterns
     *
     * A constraint without a finite bound gets no row.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::allocate()
    {
      if (verify() != 0)
      {
        Log::error() << "OptimalPowerFlow::SystemModel: invalid model data\n";
        return 1;
      }

      bool has_infinite_bus = false;
      for (const auto& [number, bus] : buses_)
      {
        has_infinite_bus = has_infinite_bus || bus->infinite();
      }
      if (!has_infinite_bus && !buses_.empty())
      {
        buses_.begin()->second->setReference();
      }

      IdxT n = 0;
      IdxT m = 0;
      for (auto& component : components_)
      {
        auto& variables = component->variableIndices();
        for (IdxT j = 0; j < component->sizeInternal(); ++j)
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

      RealT* g_lower = g_lower_.getData();
      RealT* g_upper = g_upper_.getData();
      for (const auto& component : components_)
      {
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
      g_lower_.setDataUpdated();
      g_upper_.setDataUpdated();

      if (initialize() != 0)
      {
        Log::error() << "OptimalPowerFlow::SystemModel: initialization from the state failed\n";
        return 1;
      }

      return allocateDerivatives();
    }

    /**
     * @brief Starting point and variable bounds from the state
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::initialize()
    {
      ScalarT* x       = x_.getData();
      RealT*   x_lower = x_lower_.getData();
      RealT*   x_upper = x_upper_.getData();

      int ret = 0;
      for (auto& component : components_)
      {
        ret += component->initialize(state_, x, x_lower, x_upper);
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
      x_lower_.setDataUpdated();
      x_upper_.setDataUpdated();

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
     * @brief Device currents follow from the power into each terminal bus
     */
    template <typename scalar_type, typename index_type>
    Model::StateData SystemModel<scalar_type, index_type>::solutionState() const
    {
      Model::StateData state = state_;
      const ScalarT*   x     = x_.getData();

      for (const auto& [number, bus] : buses_)
      {
        auto& voltage = state.buses[Model::busKey(number)].values;
        voltage["vr"] = static_cast<RealT>(x[bus->variableIndices()[0]]);
        voltage["vi"] = static_cast<RealT>(x[bus->variableIndices()[1]]);
      }

      for (const auto& component : components_)
      {
        const auto&          terminals = component->terminals();
        std::vector<ScalarT> power(ComponentT::TERMINAL_SIZE * terminals.size());
        component->evaluateTerminalPower(x, power.data());

        for (IdxT t = 0; t < terminals.size(); ++t)
        {
          Model::setTerminalCurrent(state,
                                    component->id(),
                                    terminals[t],
                                    t,
                                    terminals.size(),
                                    static_cast<RealT>(power[ComponentT::TERMINAL_SIZE * t]),
                                    static_cast<RealT>(power[ComponentT::TERMINAL_SIZE * t + 1]));
        }
      }

      return state;
    }

    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::verify() const
    {
      int ret = data_errors_;
      for (const auto& component : components_)
      {
        ret += component->verify();
      }
      return ret;
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

        const IdxT variable   = component.sizeInternal() + ComponentT::TERMINAL_SIZE * t;
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
