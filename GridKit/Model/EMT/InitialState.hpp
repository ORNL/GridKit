#pragma once

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <GridKit/Model/EMT/Signal/Signal.hpp>

namespace GridKit::EMT
{
  template <typename ScalarT, typename IdxT>
  class Component;

  /// Operating-point values exchanged before any component initializes its state.
  template <typename ScalarT, typename IdxT>
  class InitialState
  {
  public:
    using SignalT    = Signal<ScalarT, IdxT>;
    using RealT      = typename SignalT::RealT;
    using ComponentT = Component<ScalarT, IdxT>;
    using Values     = std::map<std::string, RealT>;

    struct Ports
    {
      std::vector<SignalT*>                inputs;
      std::multimap<std::string, SignalT*> outputs;
      std::vector<SignalT*>                targets;
    };

    explicit InitialState(RealT omega = RealT{0})
      : omega_(omega)
    {
      if (!std::isfinite(omega) || omega < RealT{0})
        throw std::invalid_argument("Initial angular frequency must be finite and nonnegative");
    }

    RealT omega() const
    {
      return omega_;
    }

    void add(ComponentT& component, std::string path, Values values, const Ports& ports)
    {
      try
      {
        component.validateInitialState(values);
      }
      catch (const std::exception& error)
      {
        throw std::invalid_argument(path + ": " + error.what());
      }
      entries_.emplace(&component, Entry{std::move(path), std::move(values), ports});
      for (const auto& [name, signal] : ports.outputs)
      {
        if (!owners_.emplace(signal, Output{&component, name}).second)
          throw std::invalid_argument("Multiple initializers for signal: " + signal->id());
      }
    }

    const Values& outputs(ComponentT& component) const
    {
      return entries_.at(&component).values;
    }

    ComponentT* owner(const SignalT& signal) const
    {
      const auto entry = owners_.find(&signal);
      return entry == owners_.end() ? nullptr : entry->second.component;
    }

    /// Publish a value computed without changing the producer's state.
    void provide(const SignalT& signal, RealT value)
    {
      if (!std::isfinite(value))
        throw std::invalid_argument("Nonfinite initial signal: " + signal.id());
      if (signal.constant())
      {
        check(value, static_cast<RealT>(signal.read()), "constant " + signal.id());
        return;
      }
      const auto target = owners_.find(&signal);
      if (target != owners_.end())
      {
        auto&       entry      = entries_.at(target->second.component);
        const auto& name       = target->second.name;
        const auto  prescribed = entry.values.emplace(name, value).first;
        check(value, prescribed->second, entry.path + "." + name);
        const auto [first, last] = entry.ports.outputs.equal_range(name);
        for (auto port = first; port != last; ++port)
        {
          auto planned = planned_.emplace(port->second, value).first;
          check(value, planned->second, entry.path + "." + name);
        }
      }
      auto planned = planned_.emplace(&signal, value).first;
      check(value, planned->second, signal.id());
    }

    RealT value(const SignalT& signal) const
    {
      if (signal.constant())
        return static_cast<RealT>(signal.read());
      const auto entry = planned_.find(&signal);
      if (entry == planned_.end())
        throw std::invalid_argument("Initial signal value has not been provided: " + signal.id());
      return entry->second;
    }

    /// Request an operating point from the output's owner; constants are constraints.
    void require(const SignalT& signal, RealT value, ComponentT& requester)
    {
      const auto& source = entries_.at(&requester).path;
      if (signal.constant())
      {
        check(value, static_cast<RealT>(signal.read()), source + ": constant " + signal.id());
        return;
      }
      const auto target = owners_.find(&signal);
      if (target == owners_.end())
        throw std::invalid_argument(source + ": signal has no initializable output: " + signal.id());
      auto& entry      = entries_.at(target->second.component);
      auto  prescribed = entry.values.emplace(target->second.name, value).first;
      check(value, prescribed->second, source + ": " + entry.path + "." + target->second.name);
      provide(signal, value);
    }

    int initialize()
    {
      std::map<IdxT, ComponentT*> variables;
      std::vector<ComponentT*>    components;
      for (const auto& [component, entry] : entries_)
      {
        components.push_back(component);
        for (auto index : component->getVariableIndices())
          variables.emplace(index, component);
      }
      std::sort(components.begin(), components.end(), [&](auto* lhs, auto* rhs)
                { return entries_.at(lhs).path < entries_.at(rhs).path; });

      using Dependencies = std::map<ComponentT*, std::set<ComponentT*>>;
      Dependencies preparation, state;
      for (auto* component : components)
      {
        const auto& entry = entries_.at(component);
        for (const auto* signal : entry.ports.inputs)
        {
          if (!signal || signal->constant())
            continue;
          std::set<ComponentT*> producers;
          if (auto* producer = owner(*signal))
            producers.insert(producer);
          else
          {
            typename SignalT::GradientT gradient;
            signal->appendGradient(gradient);
            for (const auto& [column, coefficient] : gradient)
            {
              const auto producer = variables.find(column);
              if (producer == variables.end())
                throw std::invalid_argument(entry.path + ": initialization input has no state owner");
              producers.insert(producer->second);
            }
          }
          producers.erase(component);
          state[component].insert(producers.begin(), producers.end());
          if (std::find(entry.ports.targets.begin(), entry.ports.targets.end(), signal) == entry.ports.targets.end())
            preparation[component].insert(producers.begin(), producers.end());
        }
        for (const auto* signal : entry.ports.targets)
        {
          if (signal->constant())
            continue;
          if (auto* producer = owner(*signal))
            preparation[producer].insert(component);
        }
      }

      auto ordered = [&](const Dependencies& dependencies)
      {
        std::vector<ComponentT*>   order;
        std::map<ComponentT*, int> visited;
        auto                       visit = [&](auto&& self, ComponentT* component) -> void
        {
          auto& mark = visited[component];
          if (mark == 1)
            throw std::invalid_argument("Cyclic initialization dependency at " + entries_.at(component).path);
          if (mark == 2)
            return;
          mark = 1;
          if (const auto entry = dependencies.find(component); entry != dependencies.end())
            for (auto* dependency : entry->second)
              self(self, dependency);
          mark = 2;
          order.push_back(component);
        };
        for (auto* component : components)
          visit(visit, component);
        return order;
      };
      const auto preparation_order = ordered(preparation);
      const auto state_order       = ordered(state);

      // All operating-point requirements are reconciled before state mutation.
      for (auto* component : preparation_order)
      {
        try
        {
          component->prepareInitialization(*this);
        }
        catch (const std::exception& error)
        {
          throw std::invalid_argument(entries_.at(component).path + ": " + error.what());
        }
      }
      for (auto* component : state_order)
        component->validateInitialState(outputs(*component));
      auto initialize = [&](ComponentT* component)
      {
        try
        {
          return component->initializeState(outputs(*component), omega_);
        }
        catch (const std::exception& error)
        {
          throw std::invalid_argument(entries_.at(component).path + ": " + error.what());
        }
      };
      for (auto* component : state_order)
        if (component->size() != 0)
          if (const int status = initialize(component); status != 0)
            return status;
      // Expressions own no state and initialize after their state producers.
      for (auto* component : state_order)
        if (component->size() == 0)
          if (const int status = initialize(component); status != 0)
            return status;
      return 0;
    }

  private:
    static void check(RealT value, RealT prescribed, const std::string& target)
    {
      if (!std::isfinite(value) || std::abs(value - prescribed) > RealT{1e-10} * (RealT{1} + std::abs(prescribed)))
        throw std::invalid_argument("Initial requirement conflicts with " + target);
    }

    struct Entry
    {
      std::string path;
      Values      values;
      Ports       ports;
    };

    struct Output
    {
      ComponentT* component;
      std::string name;
    };

    RealT                            omega_;
    std::map<ComponentT*, Entry>     entries_;
    std::map<const SignalT*, Output> owners_;
    std::map<const SignalT*, RealT>  planned_;
  };
} // namespace GridKit::EMT
