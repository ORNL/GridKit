#include "System.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include <GridKit/MemoryUtilities/MemoryUtils.hpp>
#include <GridKit/Model/OPF/Branch/Branch.hpp>
#include <GridKit/Model/OPF/Bus/Bus.hpp>
#include <GridKit/Model/OPF/Generator/Generator.hpp>
#include <GridKit/Model/OPF/Load/Load.hpp>
#include <GridKit/Model/OPF/Shunt/Shunt.hpp>

namespace GridKit
{
  namespace OPF
  {
    namespace
    {
      template <typename IdxT>
      bool isNegative(IdxT value)
      {
        if constexpr (std::is_signed_v<IdxT>)
        {
          return value < 0;
        }
        return false;
      }

      template <typename IdxT>
      std::string busStateId(IdxT number)
      {
        return "bus_id_" + std::to_string(number);
      }

      [[noreturn]] void throwValidationError(const std::string& context,
                                             const std::string& message)
      {
        throw std::invalid_argument("OPF " + context + ": " + message);
      }

      template <typename RealT>
      void requireFinite(RealT              value,
                         const std::string& context)
      {
        if (!std::isfinite(value))
        {
          throwValidationError(context, "expected a finite value");
        }
      }

      template <typename RealT>
      void requirePositive(RealT              value,
                           const std::string& context)
      {
        requireFinite(value, context);
        if (!(value > RealT{0}))
        {
          throwValidationError(context, "expected a positive value");
        }
      }

      template <typename RealT>
      void validateLimits(const std::optional<RealT>& lower,
                          const std::optional<RealT>& upper,
                          const std::string&          lower_context,
                          const std::string&          upper_context)
      {
        if (lower)
        {
          requireFinite(*lower, lower_context);
        }
        if (upper)
        {
          requireFinite(*upper, upper_context);
        }
        if (lower && upper && *lower > *upper)
        {
          throwValidationError(lower_context,
                               "lower limit exceeds " + upper_context);
        }
      }

      template <typename IdxT>
      void checkedAdd(IdxT&              total,
                      IdxT               increment,
                      const std::string& context)
      {
        if (isNegative(total) || isNegative(increment)
            || total > std::numeric_limits<IdxT>::max() - increment)
        {
          throw std::overflow_error("OPF " + context
                                    + " exceeds the supported index range");
        }
        total += increment;
      }

      template <typename IdxT>
      std::size_t checkedSize(IdxT value, const std::string& context)
      {
        if (isNegative(value))
        {
          throw std::logic_error("OPF " + context + " is negative");
        }
        const auto size = static_cast<std::size_t>(value);
        if (static_cast<IdxT>(size) != value)
        {
          throw std::overflow_error("OPF " + context
                                    + " exceeds the host size range");
        }
        return size;
      }

      template <typename ValueT>
      bool allFinite(const ValueT* values, std::size_t count)
      {
        for (std::size_t i = 0; i < count; ++i)
        {
          if (!std::isfinite(values[i]))
          {
            return false;
          }
        }
        return true;
      }
    } // namespace

    template <class ScalarT, typename IdxT>
    System<ScalarT, IdxT>::System(const SystemDataT&      system_data,
                                  const Model::StateData& state)
      : system_data_(system_data),
        input_state_(state)
    {
    }

    template <class ScalarT, typename IdxT>
    void System<ScalarT, IdxT>::validate() const
    {
      requirePositive(system_data_.params.freq_base, "params.freq_base");
      requirePositive(system_data_.params.va_base, "params.va_base");

      std::set<std::string>             bus_ids;
      std::set<IdxT>                    bus_numbers;
      std::map<IdxT, std::vector<IdxT>> adjacency;
      IdxT                              slack_number{};
      std::size_t                       slack_count = 0;

      for (std::size_t i = 0; i < system_data_.buses.size(); ++i)
      {
        const auto& bus     = system_data_.buses[i];
        const auto  context = "buses[" + std::to_string(i) + "] (\""
                             + bus.id + "\")";

        if (bus.id.empty())
        {
          throwValidationError(context + ".id", "must not be empty");
        }
        if (!bus_ids.insert(bus.id).second)
        {
          throwValidationError(context + ".id", "duplicate bus id");
        }
        if (isNegative(bus.number))
        {
          throwValidationError(context + ".params.number",
                               "expected a non-negative bus number");
        }
        if (!bus_numbers.insert(bus.number).second)
        {
          throwValidationError(context + ".params.number",
                               "duplicate bus number");
        }
        if (bus.bus_class != BusClass::BUS
            && bus.bus_class != BusClass::SLACK)
        {
          throwValidationError(context + ".class", "unsupported bus class");
        }

        requirePositive(bus.kv, context + ".params.kv");
        validateLimits(bus.vmin,
                       bus.vmax,
                       context + ".params.vmin",
                       context + ".params.vmax");
        if (bus.vmin && *bus.vmin < RealT{0})
        {
          throwValidationError(context + ".params.vmin",
                               "voltage magnitude limit must be non-negative");
        }
        if (bus.vmax && !(*bus.vmax > RealT{0}))
        {
          throwValidationError(context + ".params.vmax",
                               "voltage magnitude limit must be positive");
        }

        const auto state_id = busStateId(bus.number);
        const auto state_it = input_state_.buses.find(state_id);
        if (state_it == input_state_.buses.end())
        {
          throwValidationError("state.buses." + state_id,
                               "required bus state is missing");
        }
        const auto& state = state_it->second;
        if (!state.vr || !state.vi)
        {
          throwValidationError("state.buses." + state_id,
                               "both vr and vi are required");
        }
        requireFinite(*state.vr, "state.buses." + state_id + ".vr");
        requireFinite(*state.vi, "state.buses." + state_id + ".vi");
        const RealT vm = std::hypot(static_cast<RealT>(*state.vr),
                                    static_cast<RealT>(*state.vi));
        if (!(vm > RealT{0}) || !std::isfinite(vm))
        {
          throwValidationError("state.buses." + state_id,
                               "voltage magnitude must be finite and nonzero");
        }

        adjacency.emplace(bus.number, std::vector<IdxT>{});
        if (bus.bus_class == BusClass::SLACK)
        {
          slack_number = bus.number;
          ++slack_count;
        }
      }

      if (slack_count != 1)
      {
        throwValidationError("buses",
                             "expected exactly one Slack bus; found "
                                 + std::to_string(slack_count));
      }

      std::set<std::string> device_ids;
      for (std::size_t i = 0; i < system_data_.devices.size(); ++i)
      {
        const auto&       device = system_data_.devices[i];
        const std::string id     = std::visit(
            [](const auto& data)
            { return data.id; },
            device);
        const auto context = "devices[" + std::to_string(i) + "] (\""
                             + id + "\")";

        if (id.empty())
        {
          throwValidationError(context + ".id", "must not be empty");
        }
        if (!device_ids.insert(id).second)
        {
          throwValidationError(context + ".id", "duplicate device id");
        }

        const auto requireBus = [&](IdxT number, const std::string& field)
        {
          if (isNegative(number) || !bus_numbers.contains(number))
          {
            throwValidationError(context + ".buses." + field,
                                 "references unknown bus "
                                     + std::to_string(number));
          }
        };

        const auto state_it = input_state_.devices.find(id);
        std::visit(
            [&](const auto& data)
            {
              using DataT = std::decay_t<decltype(data)>;
              if constexpr (std::is_same_v<DataT, typename SystemDataT::BranchDataT>)
              {
                requireBus(data.from, "from");
                requireBus(data.to, "to");
                if (data.from == data.to)
                {
                  throwValidationError(context + ".buses",
                                       "branch terminals must be distinct");
                }
                requireFinite(data.R, context + ".params.R");
                requireFinite(data.X, context + ".params.X");
                requireFinite(data.G, context + ".params.G");
                requireFinite(data.B, context + ".params.B");
                const RealT impedance_squared = data.R * data.R + data.X * data.X;
                if (!std::isfinite(impedance_squared)
                    || !(impedance_squared > RealT{0}))
                {
                  throwValidationError(context + ".params",
                                       "R squared plus X squared must be finite and positive");
                }
                if (data.smax)
                {
                  requireFinite(*data.smax, context + ".params.smax");
                  if (*data.smax < RealT{0})
                  {
                    throwValidationError(context + ".params.smax",
                                         "must be non-negative");
                  }
                  if (!std::isfinite(*data.smax * *data.smax))
                  {
                    throwValidationError(context + ".params.smax",
                                         "squared rating must be finite");
                  }
                }

                RealT tap   = RealT{1};
                RealT phase = RealT{0};
                bool  open  = false;
                if (state_it != input_state_.devices.end())
                {
                  tap   = state_it->second.tap.value_or(tap);
                  phase = state_it->second.phase.value_or(phase);
                  open  = state_it->second.open.value_or(open);
                }
                requirePositive(tap, "state.devices." + id + ".tap");
                requireFinite(phase, "state.devices." + id + ".phase");

                const RealT                inverse_tap         = RealT{1} / tap;
                const RealT                inverse_tap_squared = inverse_tap * inverse_tap;
                const RealT                series_G            = data.R / impedance_squared;
                const RealT                series_B            = -data.X / impedance_squared;
                const RealT                cos_phase           = std::cos(phase);
                const RealT                sin_phase           = std::sin(phase);
                const std::array<RealT, 8> admittance{
                    -(series_G + RealT{0.5} * data.G) * inverse_tap_squared,
                    -(series_B + RealT{0.5} * data.B) * inverse_tap_squared,
                    (series_G * cos_phase - series_B * sin_phase) * inverse_tap,
                    (series_B * cos_phase + series_G * sin_phase) * inverse_tap,
                    (series_G * cos_phase + series_B * sin_phase) * inverse_tap,
                    (series_B * cos_phase - series_G * sin_phase) * inverse_tap,
                    -(series_G + RealT{0.5} * data.G),
                    -(series_B + RealT{0.5} * data.B)};
                if (!std::isfinite(inverse_tap)
                    || !std::isfinite(inverse_tap_squared)
                    || !std::all_of(admittance.begin(),
                                    admittance.end(),
                                    [](RealT value)
                                    { return std::isfinite(value); }))
                {
                  throwValidationError("state.devices." + id,
                                       "branch operating inputs produce non-finite admittance");
                }
                if (!open)
                {
                  adjacency.at(data.from).push_back(data.to);
                  adjacency.at(data.to).push_back(data.from);
                }
              }
              else if constexpr (std::is_same_v<DataT, typename SystemDataT::GeneratorDataT>)
              {
                requireBus(data.bus, "bus");
                validateLimits(data.pmin,
                               data.pmax,
                               context + ".params.pmin",
                               context + ".params.pmax");
                validateLimits(data.qmin,
                               data.qmax,
                               context + ".params.qmin",
                               context + ".params.qmax");
                requirePositive(data.mva, context + ".params.mva");
                requireFinite(data.c0, context + ".params.c0");
                requireFinite(data.c1, context + ".params.c1");
                requireFinite(data.c2, context + ".params.c2");
                if (state_it != input_state_.devices.end())
                {
                  if (state_it->second.p)
                  {
                    requireFinite(*state_it->second.p,
                                  "state.devices." + id + ".p");
                  }
                  if (state_it->second.q)
                  {
                    requireFinite(*state_it->second.q,
                                  "state.devices." + id + ".q");
                  }
                }
              }
              else if constexpr (std::is_same_v<DataT, typename SystemDataT::LoadDataT>)
              {
                requireBus(data.bus, "bus");
                validateLimits(data.pmin,
                               data.pmax,
                               context + ".params.pmin",
                               context + ".params.pmax");
                validateLimits(data.qmin,
                               data.qmax,
                               context + ".params.qmin",
                               context + ".params.qmax");
                if (state_it == input_state_.devices.end()
                    || !state_it->second.p || !state_it->second.q)
                {
                  throwValidationError("state.devices." + id,
                                       "load requires finite p and q values");
                }
                const RealT p = *state_it->second.p;
                const RealT q = *state_it->second.q;
                requireFinite(p, "state.devices." + id + ".p");
                requireFinite(q, "state.devices." + id + ".q");
                if ((data.pmin && p < *data.pmin)
                    || (data.pmax && p > *data.pmax))
                {
                  throwValidationError("state.devices." + id + ".p",
                                       "load active power is outside its limits");
                }
                if ((data.qmin && q < *data.qmin)
                    || (data.qmax && q > *data.qmax))
                {
                  throwValidationError("state.devices." + id + ".q",
                                       "load reactive power is outside its limits");
                }
              }
              else if constexpr (std::is_same_v<DataT, typename SystemDataT::ShuntDataT>)
              {
                requireBus(data.bus, "bus");
                requireFinite(data.G, context + ".params.G");
                requireFinite(data.B, context + ".params.B");
              }
            },
            device);
      }

      std::set<IdxT>    visited;
      std::vector<IdxT> pending{slack_number};
      while (!pending.empty())
      {
        const IdxT current = pending.back();
        pending.pop_back();
        if (!visited.insert(current).second)
        {
          continue;
        }
        for (const IdxT neighbor : adjacency.at(current))
        {
          if (!visited.contains(neighbor))
          {
            pending.push_back(neighbor);
          }
        }
      }
      if (visited.size() != system_data_.buses.size())
      {
        throwValidationError("topology",
                             "closed branches do not connect every bus to the Slack bus");
      }
    }

    template <class ScalarT, typename IdxT>
    int System<ScalarT, IdxT>::allocate()
    {
      if (allocated_)
      {
        return 0;
      }

      validate();

      using BusT       = Bus<ScalarT, IdxT>;
      using GeneratorT = Generator<ScalarT, IdxT>;
      using LoadT      = Load<ScalarT, IdxT>;
      using ShuntT     = Shunt<ScalarT, IdxT>;
      using BranchT    = Branch<ScalarT, IdxT>;

      std::vector<std::unique_ptr<ComponentT>> components;
      components.reserve(system_data_.buses.size() + system_data_.devices.size());
      std::map<IdxT, BusT*> buses;
      IdxT                  variable_count{};
      IdxT                  constraint_count{};

      auto append = [&](std::unique_ptr<ComponentT> component)
      {
        const IdxT variable_size   = component->sizeInternalVariables();
        const IdxT constraint_size = component->sizeInternalConstraints();

        IdxT next_variable_count   = variable_count;
        IdxT next_constraint_count = constraint_count;
        checkedAdd(next_variable_count, variable_size, "variable count");
        checkedAdd(next_constraint_count, constraint_size, "constraint count");

        component->setVariableOffset(variable_count);
        component->setConstraintOffset(constraint_count);
        variable_count   = next_variable_count;
        constraint_count = next_constraint_count;
        components.push_back(std::move(component));
      };

      for (const auto& data : system_data_.buses)
      {
        const auto state_it = input_state_.buses.find(busStateId(data.number));
        auto       bus      = std::make_unique<BusT>(data, state_it->second);
        BusT*      bus_ptr  = bus.get();
        append(std::move(bus));
        buses.emplace(data.number, bus_ptr);
      }

      const auto deviceState = [&](const std::string& id)
      {
        Model::DeviceState state;
        const auto         entry = input_state_.devices.find(id);
        if (entry != input_state_.devices.end())
        {
          state = entry->second;
        }
        return state;
      };

      for (const auto& device : system_data_.devices)
      {
        std::visit(
            [&](const auto& data)
            {
              using DataT = std::decay_t<decltype(data)>;
              if constexpr (std::is_same_v<DataT, typename SystemDataT::BranchDataT>)
              {
                Model::DeviceState state = deviceState(data.id);
                state.open               = state.open.value_or(false);
                state.tap                = state.tap.value_or(RealT{1});
                state.phase              = state.phase.value_or(RealT{0});
                auto  branch             = std::make_unique<BranchT>(data, state);
                BusT* from               = buses.at(data.from);
                BusT* to                 = buses.at(data.to);

                branch->variables().template bindExternal<BranchExternalVariables::VMF>(
                    from->variables().template internalIndex<BusInternalVariables::VM>());
                branch->variables().template bindExternal<BranchExternalVariables::VAF>(
                    from->variables().template internalIndex<BusInternalVariables::VA>());
                branch->variables().template bindExternal<BranchExternalVariables::VMT>(
                    to->variables().template internalIndex<BusInternalVariables::VM>());
                branch->variables().template bindExternal<BranchExternalVariables::VAT>(
                    to->variables().template internalIndex<BusInternalVariables::VA>());
                branch->constraints().template bindExternal<BranchExternalConstraints::DIVPF>(
                    from->constraints().template internalIndex<BusInternalConstraints::DIVP>());
                branch->constraints().template bindExternal<BranchExternalConstraints::DIVQF>(
                    from->constraints().template internalIndex<BusInternalConstraints::DIVQ>());
                branch->constraints().template bindExternal<BranchExternalConstraints::DIVPT>(
                    to->constraints().template internalIndex<BusInternalConstraints::DIVP>());
                branch->constraints().template bindExternal<BranchExternalConstraints::DIVQT>(
                    to->constraints().template internalIndex<BusInternalConstraints::DIVQ>());
                append(std::move(branch));
              }
              else if constexpr (std::is_same_v<DataT, typename SystemDataT::GeneratorDataT>)
              {
                Model::DeviceState state = deviceState(data.id);
                state.online             = state.online.value_or(true);
                state.p                  = state.p.value_or(RealT{0});
                state.q                  = state.q.value_or(RealT{0});
                auto  generator          = std::make_unique<GeneratorT>(data, state);
                BusT* bus                = buses.at(data.bus);

                generator->variables().template bindExternal<GeneratorExternalVariables::VM>(
                    bus->variables().template internalIndex<BusInternalVariables::VM>());
                generator->variables().template bindExternal<GeneratorExternalVariables::VA>(
                    bus->variables().template internalIndex<BusInternalVariables::VA>());
                generator->constraints().template bindExternal<GeneratorExternalConstraints::DIVP>(
                    bus->constraints().template internalIndex<BusInternalConstraints::DIVP>());
                generator->constraints().template bindExternal<GeneratorExternalConstraints::DIVQ>(
                    bus->constraints().template internalIndex<BusInternalConstraints::DIVQ>());
                append(std::move(generator));
              }
              else if constexpr (std::is_same_v<DataT, typename SystemDataT::LoadDataT>)
              {
                Model::DeviceState state = deviceState(data.id);
                state.online             = state.online.value_or(true);
                auto  load               = std::make_unique<LoadT>(data, state);
                BusT* bus                = buses.at(data.bus);

                load->variables().template bindExternal<LoadExternalVariables::VM>(
                    bus->variables().template internalIndex<BusInternalVariables::VM>());
                load->variables().template bindExternal<LoadExternalVariables::VA>(
                    bus->variables().template internalIndex<BusInternalVariables::VA>());
                load->constraints().template bindExternal<LoadExternalConstraints::DIVP>(
                    bus->constraints().template internalIndex<BusInternalConstraints::DIVP>());
                load->constraints().template bindExternal<LoadExternalConstraints::DIVQ>(
                    bus->constraints().template internalIndex<BusInternalConstraints::DIVQ>());
                append(std::move(load));
              }
              else if constexpr (std::is_same_v<DataT, typename SystemDataT::ShuntDataT>)
              {
                Model::DeviceState state = deviceState(data.id);
                state.online             = state.online.value_or(true);
                auto  shunt              = std::make_unique<ShuntT>(data, state);
                BusT* bus                = buses.at(data.bus);

                shunt->variables().template bindExternal<ShuntExternalVariables::VM>(
                    bus->variables().template internalIndex<BusInternalVariables::VM>());
                shunt->variables().template bindExternal<ShuntExternalVariables::VA>(
                    bus->variables().template internalIndex<BusInternalVariables::VA>());
                shunt->constraints().template bindExternal<ShuntExternalConstraints::DIVP>(
                    bus->constraints().template internalIndex<BusInternalConstraints::DIVP>());
                shunt->constraints().template bindExternal<ShuntExternalConstraints::DIVQ>(
                    bus->constraints().template internalIndex<BusInternalConstraints::DIVQ>());
                append(std::move(shunt));
              }
            },
            device);
      }

      int status = variables_.resize(variable_count);
      if (status != 0)
      {
        return status;
      }
      status = constraints_.resize(constraint_count);
      if (status != 0)
      {
        return status;
      }
      status = objective_gradient_.resize(variable_count);
      if (status != 0)
      {
        return status;
      }
      status = variable_lower_.resize(variable_count);
      if (status != 0)
      {
        return status;
      }
      status = variable_upper_.resize(variable_count);
      if (status != 0)
      {
        return status;
      }
      status = constraint_lower_.resize(constraint_count);
      if (status != 0)
      {
        return status;
      }
      status = constraint_upper_.resize(constraint_count);
      if (status != 0)
      {
        return status;
      }

      const RealT infinity  = std::numeric_limits<RealT>::infinity();
      status                = variables_.setToZero();
      status               |= constraints_.setToZero();
      status               |= objective_gradient_.setToZero();
      status               |= variable_lower_.setToConst(-infinity);
      status               |= variable_upper_.setToConst(infinity);
      status               |= constraint_lower_.setToZero();
      status               |= constraint_upper_.setToZero();
      if (status != 0)
      {
        return status;
      }

      RealT* variable_lower   = variable_lower_.getData();
      RealT* variable_upper   = variable_upper_.getData();
      RealT* constraint_lower = constraint_lower_.getData();
      RealT* constraint_upper = constraint_upper_.getData();
      for (const auto& component : components)
      {
        status = component->setVariableBounds(variable_lower, variable_upper);
        if (status != 0)
        {
          return status;
        }
        status = component->setConstraintBounds(constraint_lower, constraint_upper);
        if (status != 0)
        {
          return status;
        }
      }
      variable_lower_.setDataUpdated();
      variable_upper_.setDataUpdated();
      constraint_lower_.setDataUpdated();
      constraint_upper_.setDataUpdated();

      using JacobianEntry = typename ComponentT::JacobianEntry;
      std::vector<JacobianEntry>                       entries;
      std::vector<std::pair<std::size_t, std::size_t>> component_spans;
      component_spans.reserve(components.size());
      for (const auto& component : components)
      {
        const std::size_t begin = entries.size();
        component->addJacobianSparsity(entries);
        const std::size_t count = entries.size() - begin;
        if (count != checkedSize(component->nnz(), "component Jacobian count"))
        {
          throw std::logic_error(
              "OPF component Jacobian sparsity count does not match nnz()");
        }
        component_spans.emplace_back(begin, count);
      }

      for (const auto& [row, column] : entries)
      {
        if (isNegative(row) || row >= constraint_count
            || isNegative(column) || column >= variable_count)
        {
          throw std::logic_error("OPF component returned an invalid Jacobian coordinate");
        }
      }

      std::vector<std::size_t> order(entries.size());
      std::iota(order.begin(), order.end(), std::size_t{0});
      std::sort(order.begin(),
                order.end(),
                [&](std::size_t left, std::size_t right)
                {
                  if (entries[left].first != entries[right].first)
                  {
                    return entries[left].first < entries[right].first;
                  }
                  if (entries[left].second != entries[right].second)
                  {
                    return entries[left].second < entries[right].second;
                  }
                  return left < right;
                });

      std::vector<JacobianEntry> unique_entries;
      std::vector<IdxT>          raw_to_slot(entries.size());
      unique_entries.reserve(entries.size());
      for (const std::size_t raw_index : order)
      {
        if (unique_entries.empty() || unique_entries.back() != entries[raw_index])
        {
          unique_entries.push_back(entries[raw_index]);
        }
        const std::size_t slot = unique_entries.size() - 1;
        if (slot > static_cast<std::size_t>(std::numeric_limits<IdxT>::max()))
        {
          throw std::overflow_error("OPF Jacobian exceeds the supported index range");
        }
        raw_to_slot[raw_index] = static_cast<IdxT>(slot);
      }

      if (unique_entries.size()
          > static_cast<std::size_t>(std::numeric_limits<IdxT>::max()))
      {
        throw std::overflow_error("OPF Jacobian exceeds the supported index range");
      }
      const IdxT unique_nnz = static_cast<IdxT>(unique_entries.size());
      auto       jacobian   = std::make_unique<CsrMatrixT>(constraint_count,
                                                   variable_count,
                                                   unique_nnz);
      status                = jacobian->allocateMatrixData(memory::HOST);
      if (status != 0)
      {
        return status;
      }

      IdxT* row_data = jacobian->getRowData();
      IdxT* col_data = jacobian->getColData();
      for (const auto& [row, column] : unique_entries)
      {
        static_cast<void>(column);
        ++row_data[static_cast<std::size_t>(row) + 1];
      }
      for (IdxT row = 0; row < constraint_count; ++row)
      {
        row_data[static_cast<std::size_t>(row) + 1] += row_data[row];
      }
      for (std::size_t slot = 0; slot < unique_entries.size(); ++slot)
      {
        col_data[slot] = unique_entries[slot].second;
      }
      jacobian->setUpdated(memory::HOST);

      for (std::size_t i = 0; i < components.size(); ++i)
      {
        const auto [begin, count] = component_spans[i];
        std::vector<IdxT> slots(raw_to_slot.begin() + static_cast<std::ptrdiff_t>(begin),
                                raw_to_slot.begin()
                                    + static_cast<std::ptrdiff_t>(begin + count));
        status = components[i]->setJacobianSlots(slots);
        if (status != 0)
        {
          throw std::logic_error("OPF component rejected its Jacobian slot mapping");
        }
      }

      components_.swap(components);
      jacobian_         = std::move(jacobian);
      size_variables_   = variable_count;
      size_constraints_ = constraint_count;
      nnz_              = unique_nnz;
      objective_        = RealT{0};
      initialized_      = false;
      allocated_        = true;
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int System<ScalarT, IdxT>::initialize()
    {
      if (!allocated_)
      {
        throw std::logic_error("OPF System must be allocated before initialization");
      }

      initialized_  = false;
      int status    = variables_.setToZero();
      status       |= constraints_.setToZero();
      status       |= objective_gradient_.setToZero();
      if (status != 0)
      {
        return status;
      }

      ScalarT* values = variables_.getData();
      for (const auto& component : components_)
      {
        status = component->initialize(values);
        if (status != 0)
        {
          return status;
        }
      }

      RealT* jacobian_values = jacobian_->getValues();
      if (nnz_ != IdxT{})
      {
        std::fill_n(jacobian_values, checkedSize(nnz_, "Jacobian size"), RealT{0});
      }
      status      = variables_.setDataUpdated();
      status     |= constraints_.setDataUpdated();
      status     |= objective_gradient_.setDataUpdated();
      status     |= jacobian_->setUpdated(memory::HOST);
      objective_  = RealT{0};
      if (status == 0)
      {
        initialized_ = true;
      }
      return status;
    }

    template <class ScalarT, typename IdxT>
    int System<ScalarT, IdxT>::evaluateResidual()
    {
      if (!allocated_ || !initialized_)
      {
        throw std::logic_error("OPF System must be initialized before residual evaluation");
      }

      int status = constraints_.setToZero();
      if (status != 0)
      {
        return status;
      }
      const ScalarT* variables   = variables_.getData();
      ScalarT*       constraints = constraints_.getData();
      for (const auto& component : components_)
      {
        status = component->addConstraints(variables, constraints);
        if (status != 0)
        {
          return status;
        }
      }
      if (!allFinite(constraints,
                     checkedSize(size_constraints_, "constraint count")))
      {
        return 1;
      }
      return constraints_.setDataUpdated();
    }

    template <class ScalarT, typename IdxT>
    int System<ScalarT, IdxT>::evaluateJacobian()
    {
      if (!allocated_ || !initialized_)
      {
        throw std::logic_error("OPF System must be initialized before Jacobian evaluation");
      }

      int status = objective_gradient_.setToZero();
      if (status != 0)
      {
        return status;
      }
      RealT* jacobian_values = jacobian_->getValues();
      if (nnz_ != IdxT{})
      {
        std::fill_n(jacobian_values, checkedSize(nnz_, "Jacobian size"), RealT{0});
      }

      const ScalarT* variables = variables_.getData();
      ScalarT*       gradient  = objective_gradient_.getData();
      for (const auto& component : components_)
      {
        status = component->addObjectiveGradient(variables, gradient);
        if (status != 0)
        {
          return status;
        }
        status = component->addJacobian(variables, jacobian_values);
        if (status != 0)
        {
          return status;
        }
      }

      if (!allFinite(gradient,
                     checkedSize(size_variables_, "variable count"))
          || !allFinite(jacobian_values,
                        checkedSize(nnz_, "Jacobian size")))
      {
        return 1;
      }

      status  = objective_gradient_.setDataUpdated();
      status |= jacobian_->setUpdated(memory::HOST);
      return status;
    }

    template <class ScalarT, typename IdxT>
    int System<ScalarT, IdxT>::evaluateObjective()
    {
      if (!allocated_ || !initialized_)
      {
        throw std::logic_error("OPF System must be initialized before objective evaluation");
      }

      objective_               = RealT{0};
      const ScalarT* variables = variables_.getData();
      for (const auto& component : components_)
      {
        const int status = component->addObjective(variables, objective_);
        if (status != 0)
        {
          return status;
        }
      }
      if (!std::isfinite(objective_))
      {
        return 1;
      }
      return 0;
    }

    template <class ScalarT, typename IdxT>
    IdxT System<ScalarT, IdxT>::size()
    {
      return size_variables_;
    }

    template <class ScalarT, typename IdxT>
    IdxT System<ScalarT, IdxT>::sizeResidual()
    {
      return size_constraints_;
    }

    template <class ScalarT, typename IdxT>
    IdxT System<ScalarT, IdxT>::nnz()
    {
      return nnz_;
    }

    template <class ScalarT, typename IdxT>
    int System<ScalarT, IdxT>::getVariableBounds(RealVectorT& lower,
                                                 RealVectorT& upper)
    {
      if (!allocated_)
      {
        throw std::logic_error("OPF System must be allocated before requesting bounds");
      }
      int status = lower.resize(size_variables_);
      if (status != 0)
      {
        return status;
      }
      status = upper.resize(size_variables_);
      if (status != 0)
      {
        return status;
      }
      status  = lower.copyFromExternal(variable_lower_);
      status |= upper.copyFromExternal(variable_upper_);
      return status;
    }

    template <class ScalarT, typename IdxT>
    int System<ScalarT, IdxT>::getResidualBounds(RealVectorT& lower,
                                                 RealVectorT& upper)
    {
      if (!allocated_)
      {
        throw std::logic_error("OPF System must be allocated before requesting bounds");
      }
      int status = lower.resize(size_constraints_);
      if (status != 0)
      {
        return status;
      }
      status = upper.resize(size_constraints_);
      if (status != 0)
      {
        return status;
      }
      status  = lower.copyFromExternal(constraint_lower_);
      status |= upper.copyFromExternal(constraint_upper_);
      return status;
    }

    template <class ScalarT, typename IdxT>
    bool System<ScalarT, IdxT>::hasObjective() const
    {
      return true;
    }

    template <class ScalarT, typename IdxT>
    typename System<ScalarT, IdxT>::RealT System<ScalarT, IdxT>::objective() const
    {
      return objective_;
    }

    template <class ScalarT, typename IdxT>
    bool System<ScalarT, IdxT>::hasJacobian()
    {
      return jacobian_ != nullptr;
    }

    template <class ScalarT, typename IdxT>
    typename System<ScalarT, IdxT>::CsrMatrixT* System<ScalarT, IdxT>::getCsrJacobian() const
    {
      return jacobian_.get();
    }

    template <class ScalarT, typename IdxT>
    typename System<ScalarT, IdxT>::VectorT& System<ScalarT, IdxT>::y()
    {
      return variables_;
    }

    template <class ScalarT, typename IdxT>
    const typename System<ScalarT, IdxT>::VectorT& System<ScalarT, IdxT>::y() const
    {
      return variables_;
    }

    template <class ScalarT, typename IdxT>
    typename System<ScalarT, IdxT>::VectorT& System<ScalarT, IdxT>::getResidual()
    {
      return constraints_;
    }

    template <class ScalarT, typename IdxT>
    const typename System<ScalarT, IdxT>::VectorT& System<ScalarT, IdxT>::getResidual() const
    {
      return constraints_;
    }

    template <class ScalarT, typename IdxT>
    typename System<ScalarT, IdxT>::VectorT& System<ScalarT, IdxT>::getObjectiveGradient()
    {
      return objective_gradient_;
    }

    template <class ScalarT, typename IdxT>
    const typename System<ScalarT, IdxT>::VectorT& System<ScalarT, IdxT>::getObjectiveGradient() const
    {
      return objective_gradient_;
    }

    template <class ScalarT, typename IdxT>
    const typename System<ScalarT, IdxT>::SystemDataT& System<ScalarT, IdxT>::systemData() const
    {
      return system_data_;
    }

    template <class ScalarT, typename IdxT>
    const Model::StateData& System<ScalarT, IdxT>::inputState() const
    {
      return input_state_;
    }

    template <class ScalarT, typename IdxT>
    Model::StateData System<ScalarT, IdxT>::solutionState() const
    {
      if (!allocated_ || !initialized_)
      {
        throw std::logic_error("OPF System must be initialized before writing solution state");
      }

      const ScalarT*   variables = variables_.getData();
      Model::StateData state     = input_state_;
      for (const auto& component : components_)
      {
        if (component->updateSolutionState(variables, state) != 0)
        {
          throw std::runtime_error("OPF component failed to write solution state");
        }
      }
      return state;
    }

    template class System<double, int>;
    template class System<double, long int>;
    template class System<double, std::size_t>;

  } // namespace OPF
} // namespace GridKit
