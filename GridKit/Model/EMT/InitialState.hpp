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
      std::vector<SignalT*>           inputs;
      std::map<std::string, SignalT*> outputs;
      std::vector<SignalT*>           targets;
    };

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
      planned_[&signal] = value;
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
      auto        check  = [&](RealT prescribed, const std::string& target)
      {
        if (!std::isfinite(value) || std::abs(value - prescribed) > RealT{1e-10} * (RealT{1} + std::abs(prescribed)))
          throw std::invalid_argument(source + ": initial requirement conflicts with " + target);
      };
      if (signal.constant())
      {
        check(static_cast<RealT>(signal.read()), "constant " + signal.id());
        return;
      }
      const auto target = owners_.find(&signal);
      if (target == owners_.end())
        throw std::invalid_argument(source + ": signal has no initializable output: " + signal.id());
      auto& entry      = entries_.at(target->second.component);
      auto  prescribed = entry.values.emplace(target->second.name, value).first;
      check(prescribed->second, entry.path + "." + target->second.name);
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

      std::map<ComponentT*, std::set<ComponentT*>> dependencies;
      for (auto* component : components)
      {
        const auto& entry = entries_.at(component);
        for (const auto* signal : entry.ports.inputs)
        {
          if (!signal)
            continue;
          typename SignalT::GradientT gradient;
          signal->appendGradient(gradient);
          for (const auto& [column, coefficient] : gradient)
          {
            const auto producer = variables.find(column);
            if (producer == variables.end())
              throw std::invalid_argument(entry.path + ": initialization input has no state owner");
            dependencies[component].insert(producer->second);
          }
        }
        for (const auto* signal : entry.ports.targets)
        {
          if (signal->constant())
            continue;
          auto* producer = owner(*signal);
          if (!producer)
            throw std::invalid_argument(entry.path + ": signal has no initializable output: " + signal->id());
          dependencies[producer].insert(component);
        }
      }

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
        for (auto* dependency : dependencies[component])
          self(self, dependency);
        mark = 2;
        order.push_back(component);
      };
      for (auto* component : components)
        visit(visit, component);

      // All operating-point requirements are reconciled before state mutation.
      for (auto* component : order)
        component->prepareInitialization(*this);
      for (auto* component : order)
        component->validateInitialState(outputs(*component));
      auto initialize = [&](ComponentT* component)
      {
        try
        {
          return component->initializeState(outputs(*component));
        }
        catch (const std::exception& error)
        {
          throw std::invalid_argument(entries_.at(component).path + ": " + error.what());
        }
      };
      for (auto* component : order)
        if (component->size() != 0)
          if (const int status = initialize(component); status != 0)
            return status;
      // Expressions own no state; validate their prescribed values after producers initialize.
      for (auto* component : order)
        if (component->size() == 0)
          if (const int status = initialize(component); status != 0)
            return status;
      return 0;
    }

  private:
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

    std::map<ComponentT*, Entry>     entries_;
    std::map<const SignalT*, Output> owners_;
    std::map<const SignalT*, RealT>  planned_;
  };
} // namespace GridKit::EMT
