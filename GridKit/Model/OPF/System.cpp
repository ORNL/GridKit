#include "System.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <span>
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
#include <GridKit/Optimization/DerivativeStructure.hpp>

namespace GridKit
{
  namespace OPF
  {
    namespace
    {
      template <typename IdxT>
      constexpr bool isNegative(IdxT value)
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
      void requireFinite(RealT value, const std::string& context)
      {
        if (!std::isfinite(value))
        {
          throwValidationError(context, "expected a finite value");
        }
      }

      template <typename RealT>
      void requirePositive(RealT value, const std::string& context)
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
      void checkedAdd(IdxT& total, IdxT increment, const std::string& context)
      {
        if (isNegative(total) || isNegative(increment)
            || total > std::numeric_limits<IdxT>::max() - increment)
        {
          throw std::invalid_argument("OPF " + context
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
        if (values == nullptr && count != 0)
        {
          return false;
        }
        for (std::size_t entry = 0; entry < count; ++entry)
        {
          if (!std::isfinite(values[entry]))
          {
            return false;
          }
        }
        return true;
      }

      template <typename IdxT>
      struct BusBinding
      {
        std::array<IdxT, 2> variables;
        std::array<IdxT, 2> constraints;
      };
    } // namespace

    template <class ScalarT, typename IdxT>
    System<ScalarT, IdxT>::System(const SystemDataT&      system_data,
                                  const Model::StateData& state)
      : system_data_(system_data),
        input_state_(state)
    {
      static_assert(std::is_same_v<ScalarT, double>,
                    "OPF exact Enzyme derivatives currently support double only");
      static_assert(std::is_integral_v<IdxT>,
                    "OPF requires an integral index type");
    }

    template <class ScalarT, typename IdxT>
    void System<ScalarT, IdxT>::validate() const
    {
      requirePositive(system_data_.params.freq_base, "params.freq_base");
      requirePositive(system_data_.params.va_base, "params.va_base");

      if (system_data_.buses.empty())
      {
        throwValidationError("buses", "at least one bus is required");
      }

      std::set<std::string>             bus_ids;
      std::set<IdxT>                    bus_numbers;
      std::map<IdxT, std::vector<IdxT>> adjacency;

      for (std::size_t index = 0; index < system_data_.buses.size(); ++index)
      {
        const auto& bus     = system_data_.buses[index];
        const auto  context = "buses[" + std::to_string(index) + "] (\""
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
        if (!state_it->second.vr || !state_it->second.vi)
        {
          throwValidationError("state.buses." + state_id,
                               "both vr and vi are required");
        }
        requireFinite(*state_it->second.vr,
                      "state.buses." + state_id + ".vr");
        requireFinite(*state_it->second.vi,
                      "state.buses." + state_id + ".vi");
        const RealT magnitude = std::hypot(
            static_cast<RealT>(*state_it->second.vr),
            static_cast<RealT>(*state_it->second.vi));
        if (!(magnitude > RealT{0}) || !std::isfinite(magnitude))
        {
          throwValidationError("state.buses." + state_id,
                               "voltage magnitude must be finite and nonzero");
        }

        adjacency.emplace(bus.number, std::vector<IdxT>{});
      }

      std::set<std::string> device_ids;
      for (std::size_t index = 0; index < system_data_.devices.size(); ++index)
      {
        const auto& device = system_data_.devices[index];
        const auto  id     = std::visit(
            [](const auto& data)
            { return data.id; },
            device);
        const auto context = "devices[" + std::to_string(index) + "] (\""
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
              if constexpr (std::is_same_v<DataT,
                                           typename SystemDataT::BranchDataT>)
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
                if (data.R == RealT{0} && data.X == RealT{0})
                {
                  throwValidationError(context + ".params",
                                       "branch impedance must be nonzero");
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

                if (!open)
                {
                  adjacency.at(data.from).push_back(data.to);
                  adjacency.at(data.to).push_back(data.from);
                }
              }
              else if constexpr (std::is_same_v<
                                     DataT,
                                     typename SystemDataT::GeneratorDataT>)
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
              else if constexpr (std::is_same_v<DataT,
                                                typename SystemDataT::LoadDataT>)
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
              else if constexpr (std::is_same_v<DataT,
                                                typename SystemDataT::ShuntDataT>)
              {
                requireBus(data.bus, "bus");
                requireFinite(data.G, context + ".params.G");
                requireFinite(data.B, context + ".params.B");
              }
            },
            device);
      }

      std::set<IdxT>    visited;
      std::vector<IdxT> pending{system_data_.buses.front().number};
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
                             "closed branches do not connect all buses");
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
      using BranchT    = Branch<ScalarT, IdxT>;
      using GeneratorT = Generator<ScalarT, IdxT>;
      using LoadT      = Load<ScalarT, IdxT>;
      using ShuntT     = Shunt<ScalarT, IdxT>;

      std::vector<std::unique_ptr<ComponentT>> components;
      components.reserve(system_data_.buses.size() + system_data_.devices.size());

      std::map<IdxT, BusBinding<IdxT>> buses;
      IdxT                             variable_count{0};
      IdxT                             constraint_count{0};

      for (const auto& data : system_data_.buses)
      {
        const IdxT first_variable   = variable_count;
        const IdxT first_constraint = constraint_count;
        checkedAdd(variable_count, IdxT{2}, "variable count");
        checkedAdd(constraint_count, IdxT{2}, "constraint count");
        const std::array<IdxT, 2> variable_indices{{first_variable,
                                                    static_cast<IdxT>(first_variable + IdxT{1})}};
        const std::array<IdxT, 2> constraint_indices{{first_constraint,
                                                      static_cast<IdxT>(first_constraint + IdxT{1})}};

        const auto& state = input_state_.buses.at(busStateId(data.number));
        components.push_back(std::make_unique<BusT>(data,
                                                    state,
                                                    variable_indices,
                                                    constraint_indices));
        buses.emplace(data.number,
                      BusBinding<IdxT>{variable_indices, constraint_indices});
      }

      struct ThermalBound
      {
        IdxT  row;
        RealT upper;
      };

      std::vector<ThermalBound> thermal_bounds;

      const auto deviceState = [&](const std::string& id)
      {
        const auto entry = input_state_.devices.find(id);
        if (entry == input_state_.devices.end())
        {
          return Model::DeviceState{};
        }
        return entry->second;
      };

      for (const auto& device : system_data_.devices)
      {
        std::visit(
            [&](const auto& data)
            {
              using DataT      = std::decay_t<decltype(data)>;
              const auto state = deviceState(data.id);
              if constexpr (std::is_same_v<DataT,
                                           typename SystemDataT::BranchDataT>)
              {
                const auto&               from = buses.at(data.from);
                const auto&               to   = buses.at(data.to);
                const std::array<IdxT, 4> variable_indices{{from.variables[0],
                                                            from.variables[1],
                                                            to.variables[0],
                                                            to.variables[1]}};
                std::array<IdxT, 6>       constraint_indices{{from.constraints[0],
                                                              from.constraints[1],
                                                              to.constraints[0],
                                                              to.constraints[1],
                                                              IdxT{0},
                                                              IdxT{0}}};
                std::size_t               active_constraint_count = 4;
                if (data.smax)
                {
                  constraint_indices[4] = constraint_count;
                  checkedAdd(constraint_count, IdxT{1}, "constraint count");
                  constraint_indices[5] = constraint_count;
                  checkedAdd(constraint_count, IdxT{1}, "constraint count");
                  active_constraint_count = 6;
                  const RealT upper       = *data.smax * *data.smax;
                  thermal_bounds.push_back({constraint_indices[4], upper});
                  thermal_bounds.push_back({constraint_indices[5], upper});
                }
                components.push_back(std::make_unique<BranchT>(
                    data,
                    state,
                    variable_indices,
                    std::span<const IdxT>(constraint_indices.data(),
                                          active_constraint_count)));
              }
              else if constexpr (std::is_same_v<
                                     DataT,
                                     typename SystemDataT::GeneratorDataT>)
              {
                const auto& bus            = buses.at(data.bus);
                const IdxT  first_variable = variable_count;
                checkedAdd(variable_count, IdxT{2}, "variable count");
                const std::array<IdxT, 2> variable_indices{{first_variable,
                                                            static_cast<IdxT>(first_variable + IdxT{1})}};
                components.push_back(std::make_unique<GeneratorT>(
                    data,
                    state,
                    variable_indices,
                    bus.constraints,
                    bus.variables));
              }
              else if constexpr (std::is_same_v<DataT,
                                                typename SystemDataT::LoadDataT>)
              {
                const auto& bus = buses.at(data.bus);
                components.push_back(std::make_unique<LoadT>(
                    data,
                    state,
                    bus.constraints,
                    bus.variables));
              }
              else if constexpr (std::is_same_v<DataT,
                                                typename SystemDataT::ShuntDataT>)
              {
                const auto& bus = buses.at(data.bus);
                components.push_back(std::make_unique<ShuntT>(
                    data,
                    state,
                    bus.variables[0],
                    bus.constraints,
                    bus.variables));
              }
            },
            device);
      }

      int status  = variables_.resize(variable_count);
      status     |= variable_lower_bounds_.resize(variable_count);
      status     |= variable_upper_bounds_.resize(variable_count);
      status     |= objective_gradient_.resize(variable_count);
      status     |= constraints_.resize(constraint_count);
      status     |= constraint_lower_bounds_.resize(constraint_count);
      status     |= constraint_upper_bounds_.resize(constraint_count);
      if (status != 0)
      {
        return status;
      }

      const RealT infinity  = std::numeric_limits<RealT>::infinity();
      status                = variables_.setToZero(memory::HOST);
      status               |= objective_gradient_.setToZero(memory::HOST);
      status               |= constraints_.setToZero(memory::HOST);
      status               |= variable_lower_bounds_.setToConst(-infinity, memory::HOST);
      status               |= variable_upper_bounds_.setToConst(infinity, memory::HOST);
      status               |= constraint_lower_bounds_.setToZero(memory::HOST);
      status               |= constraint_upper_bounds_.setToZero(memory::HOST);
      if (status != 0)
      {
        return status;
      }

      auto* lower = variable_lower_bounds_.getData(memory::HOST);
      auto* upper = variable_upper_bounds_.getData(memory::HOST);
      for (const auto& component : components)
      {
        if (component->setVariableBounds(lower, upper) != 0)
        {
          return 1;
        }
      }

      auto* constraint_upper = constraint_upper_bounds_.getData(memory::HOST);
      for (const auto& bound : thermal_bounds)
      {
        constraint_upper[static_cast<std::size_t>(bound.row)] = bound.upper;
      }

      const auto  gauge       = buses.begin();
      const auto& gauge_state = input_state_.buses.at(busStateId(gauge->first));
      const RealT gauge_angle = std::atan2(
          static_cast<RealT>(*gauge_state.vi),
          static_cast<RealT>(*gauge_state.vr));
      const auto angle_index                       = gauge->second.variables[1];
      lower[static_cast<std::size_t>(angle_index)] = gauge_angle;
      upper[static_cast<std::size_t>(angle_index)] = gauge_angle;

      status  = variable_lower_bounds_.setDataUpdated(memory::HOST);
      status |= variable_upper_bounds_.setDataUpdated(memory::HOST);
      status |= constraint_lower_bounds_.setDataUpdated(memory::HOST);
      status |= constraint_upper_bounds_.setDataUpdated(memory::HOST);
      if (status != 0)
      {
        return status;
      }

      using SparseEntryT = Optimization::SparseEntry<IdxT>;
      std::vector<std::vector<SparseEntryT>> jacobian_entries(components.size());
      std::vector<std::vector<SparseEntryT>> hessian_entries(components.size());
      for (std::size_t component_index = 0;
           component_index < components.size();
           ++component_index)
      {
        const auto& component          = *components[component_index];
        const auto  variable_indices   = component.variableIndices();
        const auto  constraint_indices = component.constraintIndices();
        if (!Optimization::hasUniqueIndices<IdxT>(variable_indices))
        {
          throw std::logic_error("OPF component variable gather is not unique");
        }

        auto& mapped_jacobian = jacobian_entries[component_index];
        mapped_jacobian.reserve(component.jacobianEntries().size());
        for (const auto& entry : component.jacobianEntries())
        {
          if (isNegative(entry.variable)
              || static_cast<std::size_t>(entry.variable) >= variable_indices.size()
              || isNegative(entry.constraint)
              || static_cast<std::size_t>(entry.constraint) >= constraint_indices.size())
          {
            throw std::logic_error("OPF component returned an invalid Jacobian descriptor");
          }
          mapped_jacobian.push_back({constraint_indices[static_cast<std::size_t>(entry.constraint)],
                                     variable_indices[static_cast<std::size_t>(entry.variable)]});
        }

        auto& mapped_hessian = hessian_entries[component_index];
        mapped_hessian.reserve(component.hessianEntries().size());
        for (const auto& entry : component.hessianEntries())
        {
          if (isNegative(entry.row)
              || static_cast<std::size_t>(entry.row) >= variable_indices.size()
              || isNegative(entry.column)
              || static_cast<std::size_t>(entry.column) >= variable_indices.size())
          {
            throw std::logic_error("OPF component returned an invalid Hessian descriptor");
          }
          mapped_hessian.push_back(Optimization::lowerTriangle(
              variable_indices[static_cast<std::size_t>(entry.row)],
              variable_indices[static_cast<std::size_t>(entry.column)]));
        }
      }

      std::vector<typename AssemblerT::EntrySpan> jacobian_groups;
      std::vector<typename AssemblerT::EntrySpan> hessian_groups;
      jacobian_groups.reserve(components.size());
      hessian_groups.reserve(components.size());
      for (std::size_t index = 0; index < components.size(); ++index)
      {
        jacobian_groups.emplace_back(jacobian_entries[index]);
        hessian_groups.emplace_back(hessian_entries[index]);
      }

      std::vector<ContributionT> jacobian_contributions(components.size());
      std::vector<ContributionT> hessian_contributions(components.size());
      status = jacobian_assembler_.allocate(constraint_count,
                                            variable_count,
                                            jacobian_groups,
                                            jacobian_contributions);
      if (status != 0)
      {
        return status;
      }
      status = hessian_assembler_.allocate(variable_count,
                                           variable_count,
                                           hessian_groups,
                                           hessian_contributions);
      if (status != 0)
      {
        return status;
      }

      components_             = std::move(components);
      jacobian_contributions_ = std::move(jacobian_contributions);
      hessian_contributions_  = std::move(hessian_contributions);
      local_multiplier_buffers_.clear();
      local_multiplier_buffers_.reserve(components_.size());
      for (const auto& component : components_)
      {
        local_multiplier_buffers_.emplace_back(component->constraintIndices().size(),
                                               RealT{0});
      }

      variable_count_   = variable_count;
      constraint_count_ = constraint_count;
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
      int status    = variables_.setToZero(memory::HOST);
      status       |= objective_gradient_.setToZero(memory::HOST);
      status       |= constraints_.setToZero(memory::HOST);
      if (status != 0)
      {
        return status;
      }

      auto* values = variables_.getData(memory::HOST);
      for (const auto& component : components_)
      {
        if (component->initialize(values) != 0)
        {
          return 1;
        }
      }

      status      = variables_.setDataUpdated(memory::HOST);
      status     |= objective_gradient_.setDataUpdated(memory::HOST);
      status     |= constraints_.setDataUpdated(memory::HOST);
      status     |= jacobian_assembler_.clearValues();
      status     |= hessian_assembler_.clearValues();
      objective_  = RealT{0};
      if (status == 0)
      {
        initialized_ = true;
      }
      return status;
    }

    template <class ScalarT, typename IdxT>
    void System<ScalarT, IdxT>::requireInitialized(const char* operation) const
    {
      if (!allocated_ || !initialized_)
      {
        throw std::logic_error(std::string("OPF System must be initialized before ")
                               + operation);
      }
    }

    template <class ScalarT, typename IdxT>
    int System<ScalarT, IdxT>::gatherVariables()
    {
      const auto* values = variables_.getData(memory::HOST);
      if (values == nullptr)
      {
        return 1;
      }
      for (const auto& component : components_)
      {
        if (component->gatherVariables(values) != 0)
        {
          return 1;
        }
      }
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int System<ScalarT, IdxT>::evaluateObjective()
    {
      requireInitialized("objective evaluation");
      if (gatherVariables() != 0)
      {
        return 1;
      }

      objective_ = RealT{0};
      for (const auto& component : components_)
      {
        if (component->evaluateObjective() != 0)
        {
          return 1;
        }
        objective_ += component->objective();
      }
      return std::isfinite(objective_) ? 0 : 1;
    }

    template <class ScalarT, typename IdxT>
    int System<ScalarT, IdxT>::evaluateObjectiveGradient()
    {
      requireInitialized("objective-gradient evaluation");
      if (gatherVariables() != 0
          || objective_gradient_.setToZero(memory::HOST) != 0)
      {
        return 1;
      }

      auto* gradient = objective_gradient_.getData(memory::HOST);
      for (const auto& component : components_)
      {
        if (component->evaluateObjectiveGradient() != 0)
        {
          return 1;
        }
        const auto indices = component->variableIndices();
        const auto values  = component->objectiveGradientValues();
        if (indices.size() != values.size())
        {
          return 1;
        }
        for (std::size_t entry = 0; entry < values.size(); ++entry)
        {
          gradient[static_cast<std::size_t>(indices[entry])] += values[entry];
        }
      }
      if (!allFinite(gradient, checkedSize(variable_count_, "variable count")))
      {
        return 1;
      }
      return objective_gradient_.setDataUpdated(memory::HOST);
    }

    template <class ScalarT, typename IdxT>
    int System<ScalarT, IdxT>::evaluateConstraints()
    {
      requireInitialized("constraint evaluation");
      if (gatherVariables() != 0 || constraints_.setToZero(memory::HOST) != 0)
      {
        return 1;
      }

      auto* constraints = constraints_.getData(memory::HOST);
      for (const auto& component : components_)
      {
        if (component->evaluateConstraints() != 0)
        {
          return 1;
        }
        const auto indices = component->constraintIndices();
        const auto values  = component->constraintValues();
        if (indices.size() != values.size())
        {
          return 1;
        }
        for (std::size_t entry = 0; entry < values.size(); ++entry)
        {
          constraints[static_cast<std::size_t>(indices[entry])] += values[entry];
        }
      }
      if (!allFinite(constraints,
                     checkedSize(constraint_count_, "constraint count")))
      {
        return 1;
      }
      return constraints_.setDataUpdated(memory::HOST);
    }

    template <class ScalarT, typename IdxT>
    int System<ScalarT, IdxT>::evaluateJacobian()
    {
      requireInitialized("Jacobian evaluation");
      if (gatherVariables() != 0 || jacobian_assembler_.clearValues() != 0)
      {
        return 1;
      }

      for (std::size_t index = 0; index < components_.size(); ++index)
      {
        if (components_[index]->evaluateJacobian() != 0
            || jacobian_assembler_.addValues(jacobian_contributions_[index],
                                             components_[index]->jacobianValues())
                   != 0)
        {
          return 1;
        }
      }

      auto* matrix = jacobian_assembler_.matrix();
      return allFinite(matrix->getValues(memory::HOST),
                       checkedSize(matrix->getNnz(), "Jacobian size"))
                 ? 0
                 : 1;
    }

    template <class ScalarT, typename IdxT>
    int System<ScalarT, IdxT>::evaluateHessian(
        RealT        objective_factor,
        const RealT* multipliers,
        IdxT         multiplier_count)
    {
      requireInitialized("Hessian evaluation");
      if (multiplier_count != constraint_count_
          || (constraint_count_ != IdxT{0} && multipliers == nullptr)
          || !std::isfinite(objective_factor)
          || !allFinite(multipliers,
                        checkedSize(multiplier_count, "multiplier count"))
          || gatherVariables() != 0
          || hessian_assembler_.clearValues() != 0)
      {
        return 1;
      }

      for (std::size_t index = 0; index < components_.size(); ++index)
      {
        const auto constraint_indices = components_[index]->constraintIndices();
        auto&      local_multipliers  = local_multiplier_buffers_[index];
        for (std::size_t entry = 0; entry < constraint_indices.size(); ++entry)
        {
          local_multipliers[entry] =
              multipliers[static_cast<std::size_t>(constraint_indices[entry])];
        }
        if (components_[index]->evaluateHessian(objective_factor,
                                                local_multipliers)
                != 0
            || hessian_assembler_.addValues(hessian_contributions_[index],
                                            components_[index]->hessianValues())
                   != 0)
        {
          return 1;
        }
      }

      auto* matrix = hessian_assembler_.matrix();
      return allFinite(matrix->getValues(memory::HOST),
                       checkedSize(matrix->getNnz(), "Hessian size"))
                 ? 0
                 : 1;
    }

    template <class ScalarT, typename IdxT>
    IdxT System<ScalarT, IdxT>::size()
    {
      return variable_count_;
    }

    template <class ScalarT, typename IdxT>
    IdxT System<ScalarT, IdxT>::sizeConstraints()
    {
      return constraint_count_;
    }

    template <class ScalarT, typename IdxT>
    typename System<ScalarT, IdxT>::VectorT& System<ScalarT, IdxT>::variables()
    {
      return variables_;
    }

    template <class ScalarT, typename IdxT>
    const typename System<ScalarT, IdxT>::VectorT& System<ScalarT, IdxT>::variables() const
    {
      return variables_;
    }

    template <class ScalarT, typename IdxT>
    const typename System<ScalarT, IdxT>::VectorT&
    System<ScalarT, IdxT>::variableLowerBounds() const
    {
      return variable_lower_bounds_;
    }

    template <class ScalarT, typename IdxT>
    const typename System<ScalarT, IdxT>::VectorT&
    System<ScalarT, IdxT>::variableUpperBounds() const
    {
      return variable_upper_bounds_;
    }

    template <class ScalarT, typename IdxT>
    typename System<ScalarT, IdxT>::RealT System<ScalarT, IdxT>::objective() const
    {
      return objective_;
    }

    template <class ScalarT, typename IdxT>
    const typename System<ScalarT, IdxT>::VectorT&
    System<ScalarT, IdxT>::objectiveGradient() const
    {
      return objective_gradient_;
    }

    template <class ScalarT, typename IdxT>
    const typename System<ScalarT, IdxT>::VectorT& System<ScalarT, IdxT>::constraints() const
    {
      return constraints_;
    }

    template <class ScalarT, typename IdxT>
    const typename System<ScalarT, IdxT>::VectorT&
    System<ScalarT, IdxT>::constraintLowerBounds() const
    {
      return constraint_lower_bounds_;
    }

    template <class ScalarT, typename IdxT>
    const typename System<ScalarT, IdxT>::VectorT&
    System<ScalarT, IdxT>::constraintUpperBounds() const
    {
      return constraint_upper_bounds_;
    }

    template <class ScalarT, typename IdxT>
    typename System<ScalarT, IdxT>::CsrMatrixT*
    System<ScalarT, IdxT>::getCsrJacobian() const
    {
      return jacobian_assembler_.matrix();
    }

    template <class ScalarT, typename IdxT>
    typename System<ScalarT, IdxT>::CsrMatrixT*
    System<ScalarT, IdxT>::getCsrHessian() const
    {
      return hessian_assembler_.matrix();
    }

    template <class ScalarT, typename IdxT>
    bool System<ScalarT, IdxT>::hasJacobian()
    {
      return std::all_of(components_.begin(),
                         components_.end(),
                         [](const auto& component)
                         { return component->hasJacobian(); });
    }

    template <class ScalarT, typename IdxT>
    bool System<ScalarT, IdxT>::hasHessian()
    {
      return std::all_of(components_.begin(),
                         components_.end(),
                         [](const auto& component)
                         { return component->hasHessian(); });
    }

    template <class ScalarT, typename IdxT>
    const typename System<ScalarT, IdxT>::SystemDataT&
    System<ScalarT, IdxT>::systemData() const
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
      requireInitialized("solution-state output");
      const auto*      values = variables_.getData(memory::HOST);
      Model::StateData state  = input_state_;
      for (const auto& component : components_)
      {
        if (component->updateSolutionState(values, state) != 0)
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
