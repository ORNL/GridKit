#pragma once

#include <cmath>
#include <limits>
#include <string>

#include "Generator.hpp"

namespace GridKit
{
  namespace OPF
  {
    template <class ScalarT, typename IdxT>
    Generator<ScalarT, IdxT>::Generator(
        const DataT&                       data,
        const StateT&                      state,
        std::array<IdxT, VARIABLE_COUNT>   variable_indices,
        std::array<IdxT, CONSTRAINT_COUNT> constraint_indices,
        std::array<IdxT, 2>                bus_voltage_indices)
      : data_(data),
        state_(state),
        in_service_(state.online.value_or(true) ? RealT{1} : RealT{0}),
        variable_indices_(variable_indices),
        constraint_indices_(constraint_indices),
        bus_voltage_indices_(bus_voltage_indices)
    {
    }

    template <class ScalarT, typename IdxT>
    std::span<const IdxT> Generator<ScalarT, IdxT>::variableIndices() const
    {
      return variable_indices_;
    }

    template <class ScalarT, typename IdxT>
    std::span<const IdxT> Generator<ScalarT, IdxT>::constraintIndices() const
    {
      return constraint_indices_;
    }

    template <class ScalarT, typename IdxT>
    std::span<const typename Generator<ScalarT, IdxT>::JacobianEntryT>
    Generator<ScalarT, IdxT>::jacobianEntries() const
    {
      return JACOBIAN_ENTRIES;
    }

    template <class ScalarT, typename IdxT>
    std::span<const typename Generator<ScalarT, IdxT>::HessianEntryT>
    Generator<ScalarT, IdxT>::hessianEntries() const
    {
      return HESSIAN_ENTRIES;
    }

    template <class ScalarT, typename IdxT>
    int Generator<ScalarT, IdxT>::gatherVariables(
        const ScalarT* global_variables)
    {
      if (global_variables == nullptr)
      {
        return 1;
      }
      for (std::size_t local = 0; local < VARIABLE_COUNT; ++local)
      {
        variables_[local] = global_variables[variable_indices_[local]];
      }
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Generator<ScalarT, IdxT>::evaluateObjective()
    {
      const ScalarT value = evaluateObjective(variables_.data());
      if (!std::isfinite(static_cast<RealT>(value)))
      {
        return 1;
      }
      objective_ = static_cast<RealT>(value);
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Generator<ScalarT, IdxT>::evaluateConstraints()
    {
      evaluateConstraints(variables_.data(), constraint_values_.data());
      for (const ScalarT value : constraint_values_)
      {
        if (!std::isfinite(static_cast<RealT>(value)))
        {
          return 1;
        }
      }
      return 0;
    }

    template <class ScalarT, typename IdxT>
    typename Generator<ScalarT, IdxT>::RealT
    Generator<ScalarT, IdxT>::objective() const
    {
      return objective_;
    }

    template <class ScalarT, typename IdxT>
    std::span<const ScalarT>
    Generator<ScalarT, IdxT>::objectiveGradientValues() const
    {
      return objective_gradient_values_;
    }

    template <class ScalarT, typename IdxT>
    std::span<const ScalarT> Generator<ScalarT, IdxT>::constraintValues() const
    {
      return constraint_values_;
    }

    template <class ScalarT, typename IdxT>
    std::span<const typename Generator<ScalarT, IdxT>::RealT>
    Generator<ScalarT, IdxT>::jacobianValues() const
    {
      return jacobian_values_;
    }

    template <class ScalarT, typename IdxT>
    std::span<const typename Generator<ScalarT, IdxT>::RealT>
    Generator<ScalarT, IdxT>::hessianValues() const
    {
      return hessian_values_;
    }

    template <class ScalarT, typename IdxT>
    bool Generator<ScalarT, IdxT>::hasJacobian() const
    {
      return true;
    }

    template <class ScalarT, typename IdxT>
    bool Generator<ScalarT, IdxT>::hasHessian() const
    {
      return true;
    }

    template <class ScalarT, typename IdxT>
    int Generator<ScalarT, IdxT>::initialize(ScalarT* global_variables) const
    {
      if (global_variables == nullptr
          || (state_.p.has_value() && !std::isfinite(*state_.p))
          || (state_.q.has_value() && !std::isfinite(*state_.q)))
      {
        return 1;
      }

      if (in_service_ == RealT{0})
      {
        global_variables[variable_indices_[0]] = ScalarT{0};
        global_variables[variable_indices_[1]] = ScalarT{0};
      }
      else
      {
        global_variables[variable_indices_[0]] =
            static_cast<ScalarT>(state_.p.value_or(0.0));
        global_variables[variable_indices_[1]] =
            static_cast<ScalarT>(state_.q.value_or(0.0));
      }
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Generator<ScalarT, IdxT>::setVariableBounds(
        RealT* global_lower_bounds,
        RealT* global_upper_bounds) const
    {
      if (global_lower_bounds == nullptr || global_upper_bounds == nullptr)
      {
        return 1;
      }

      if (in_service_ == RealT{0})
      {
        for (const IdxT index : variable_indices_)
        {
          global_lower_bounds[index] = RealT{0};
          global_upper_bounds[index] = RealT{0};
        }
        return 0;
      }

      const RealT infinity                      = std::numeric_limits<RealT>::infinity();
      global_lower_bounds[variable_indices_[0]] = data_.pmin.value_or(-infinity);
      global_upper_bounds[variable_indices_[0]] = data_.pmax.value_or(infinity);
      global_lower_bounds[variable_indices_[1]] = data_.qmin.value_or(-infinity);
      global_upper_bounds[variable_indices_[1]] = data_.qmax.value_or(infinity);
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Generator<ScalarT, IdxT>::updateSolutionState(
        const ScalarT*    global_variables,
        Model::StateData& state) const
    {
      if (global_variables == nullptr)
      {
        return 1;
      }

      const ScalarT p = global_variables[variable_indices_[0]];
      const ScalarT q = global_variables[variable_indices_[1]];
      if (!std::isfinite(static_cast<RealT>(p))
          || !std::isfinite(static_cast<RealT>(q)))
      {
        return 1;
      }

      const auto bus_key = "bus_id_" + std::to_string(data_.bus);
      auto       bus     = state.buses.find(bus_key);
      if (bus == state.buses.end())
      {
        return 1;
      }

      ScalarT current_r{0};
      ScalarT current_i{0};
      if (in_service_ != RealT{0})
      {
        const ScalarT magnitude = global_variables[bus_voltage_indices_[0]];
        const ScalarT angle     = global_variables[bus_voltage_indices_[1]];
        const ScalarT cosine    = std::cos(angle);
        const ScalarT sine      = std::sin(angle);
        if (magnitude == ScalarT{0}
            || !std::isfinite(static_cast<RealT>(magnitude))
            || !std::isfinite(static_cast<RealT>(cosine))
            || !std::isfinite(static_cast<RealT>(sine)))
        {
          return 1;
        }
        current_r = (p * cosine + q * sine) / magnitude;
        current_i = (p * sine - q * cosine) / magnitude;
      }

      if (!std::isfinite(static_cast<RealT>(current_r))
          || !std::isfinite(static_cast<RealT>(current_i)))
      {
        return 1;
      }

      auto& device = state.devices[data_.id];
      device.p     = static_cast<double>(p);
      device.q     = static_cast<double>(q);

      auto& injection = bus->second.injections[data_.id];
      injection.ir    = static_cast<double>(current_r);
      injection.ii    = static_cast<double>(current_i);
      return 0;
    }

    template <class ScalarT, typename IdxT>
    ScalarT Generator<ScalarT, IdxT>::evaluateObjective(
        const ScalarT* local_variables) const
    {
      const ScalarT p = local_variables[0];
      return static_cast<ScalarT>(in_service_)
             * (static_cast<ScalarT>(data_.c0)
                + static_cast<ScalarT>(data_.c1) * p
                + static_cast<ScalarT>(data_.c2) * p * p);
    }

    template <class ScalarT, typename IdxT>
    void Generator<ScalarT, IdxT>::evaluateConstraints(
        const ScalarT* local_variables,
        ScalarT*       local_values) const
    {
      local_values[0] = local_variables[0];
      local_values[1] = local_variables[1];
    }
  } // namespace OPF
} // namespace GridKit
