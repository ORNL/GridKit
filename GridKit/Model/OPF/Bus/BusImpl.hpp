#pragma once

#include <cmath>
#include <limits>
#include <string>

#include "Bus.hpp"

namespace GridKit
{
  namespace OPF
  {
    template <class ScalarT, typename IdxT>
    Bus<ScalarT, IdxT>::Bus(
        const DataT&                     data,
        const StateT&                    state,
        std::array<IdxT, VARIABLE_COUNT> variable_indices,
        std::array<IdxT, 2>              balance_indices)
      : data_(data),
        state_(state),
        variable_indices_(variable_indices),
        balance_indices_(balance_indices)
    {
    }

    template <class ScalarT, typename IdxT>
    std::span<const IdxT> Bus<ScalarT, IdxT>::variableIndices() const
    {
      return variable_indices_;
    }

    template <class ScalarT, typename IdxT>
    std::span<const IdxT> Bus<ScalarT, IdxT>::constraintIndices() const
    {
      return {};
    }

    template <class ScalarT, typename IdxT>
    std::span<const IdxT> Bus<ScalarT, IdxT>::balanceIndices() const
    {
      return balance_indices_;
    }

    template <class ScalarT, typename IdxT>
    IdxT Bus<ScalarT, IdxT>::voltageMagnitudeIndex() const
    {
      return variable_indices_[0];
    }

    template <class ScalarT, typename IdxT>
    IdxT Bus<ScalarT, IdxT>::voltageAngleIndex() const
    {
      return variable_indices_[1];
    }

    template <class ScalarT, typename IdxT>
    IdxT Bus<ScalarT, IdxT>::activeBalanceIndex() const
    {
      return balance_indices_[0];
    }

    template <class ScalarT, typename IdxT>
    IdxT Bus<ScalarT, IdxT>::reactiveBalanceIndex() const
    {
      return balance_indices_[1];
    }

    template <class ScalarT, typename IdxT>
    std::span<const typename Bus<ScalarT, IdxT>::JacobianEntryT>
    Bus<ScalarT, IdxT>::jacobianEntries() const
    {
      return {};
    }

    template <class ScalarT, typename IdxT>
    std::span<const typename Bus<ScalarT, IdxT>::HessianEntryT>
    Bus<ScalarT, IdxT>::hessianEntries() const
    {
      return {};
    }

    template <class ScalarT, typename IdxT>
    int Bus<ScalarT, IdxT>::gatherVariables(const ScalarT* global_variables)
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
    int Bus<ScalarT, IdxT>::evaluateObjective()
    {
      objective_ = RealT{0};
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Bus<ScalarT, IdxT>::evaluateObjectiveGradient()
    {
      objective_gradient_values_.fill(ScalarT{0});
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Bus<ScalarT, IdxT>::evaluateConstraints()
    {
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Bus<ScalarT, IdxT>::evaluateJacobian()
    {
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Bus<ScalarT, IdxT>::evaluateHessian(
        [[maybe_unused]] RealT objective_factor,
        std::span<const RealT> local_multipliers)
    {
      return local_multipliers.empty() ? 0 : 1;
    }

    template <class ScalarT, typename IdxT>
    typename Bus<ScalarT, IdxT>::RealT Bus<ScalarT, IdxT>::objective() const
    {
      return objective_;
    }

    template <class ScalarT, typename IdxT>
    std::span<const ScalarT> Bus<ScalarT, IdxT>::objectiveGradientValues() const
    {
      return objective_gradient_values_;
    }

    template <class ScalarT, typename IdxT>
    std::span<const ScalarT> Bus<ScalarT, IdxT>::constraintValues() const
    {
      return {};
    }

    template <class ScalarT, typename IdxT>
    std::span<const typename Bus<ScalarT, IdxT>::RealT>
    Bus<ScalarT, IdxT>::jacobianValues() const
    {
      return {};
    }

    template <class ScalarT, typename IdxT>
    std::span<const typename Bus<ScalarT, IdxT>::RealT>
    Bus<ScalarT, IdxT>::hessianValues() const
    {
      return {};
    }

    template <class ScalarT, typename IdxT>
    bool Bus<ScalarT, IdxT>::hasJacobian() const
    {
      return true;
    }

    template <class ScalarT, typename IdxT>
    bool Bus<ScalarT, IdxT>::hasHessian() const
    {
      return true;
    }

    template <class ScalarT, typename IdxT>
    int Bus<ScalarT, IdxT>::initialize(ScalarT* global_variables) const
    {
      if (global_variables == nullptr || !state_.vr.has_value()
          || !state_.vi.has_value() || !std::isfinite(*state_.vr)
          || !std::isfinite(*state_.vi))
      {
        return 1;
      }

      const RealT magnitude = std::hypot(static_cast<RealT>(*state_.vr),
                                         static_cast<RealT>(*state_.vi));
      const RealT angle     = std::atan2(static_cast<RealT>(*state_.vi),
                                     static_cast<RealT>(*state_.vr));
      if (!(magnitude > RealT{0}) || !std::isfinite(magnitude)
          || !std::isfinite(angle))
      {
        return 1;
      }

      global_variables[variable_indices_[0]] = static_cast<ScalarT>(magnitude);
      global_variables[variable_indices_[1]] = static_cast<ScalarT>(angle);
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Bus<ScalarT, IdxT>::setVariableBounds(
        RealT* global_lower_bounds,
        RealT* global_upper_bounds) const
    {
      if (global_lower_bounds == nullptr || global_upper_bounds == nullptr)
      {
        return 1;
      }
      const RealT infinity                      = std::numeric_limits<RealT>::infinity();
      global_lower_bounds[variable_indices_[0]] = data_.vmin.value_or(-infinity);
      global_upper_bounds[variable_indices_[0]] = data_.vmax.value_or(infinity);
      global_lower_bounds[variable_indices_[1]] = -infinity;
      global_upper_bounds[variable_indices_[1]] = infinity;
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Bus<ScalarT, IdxT>::updateSolutionState(
        const ScalarT*    global_variables,
        Model::StateData& state) const
    {
      if (global_variables == nullptr)
      {
        return 1;
      }
      const ScalarT magnitude = global_variables[variable_indices_[0]];
      const ScalarT angle     = global_variables[variable_indices_[1]];
      const ScalarT vr        = magnitude * std::cos(angle);
      const ScalarT vi        = magnitude * std::sin(angle);
      if (!std::isfinite(static_cast<RealT>(vr))
          || !std::isfinite(static_cast<RealT>(vi)))
      {
        return 1;
      }

      const auto key = "bus_id_" + std::to_string(data_.number);
      auto       bus = state.buses.find(key);
      if (bus == state.buses.end())
      {
        return 1;
      }
      bus->second.vr = static_cast<double>(vr);
      bus->second.vi = static_cast<double>(vi);
      return 0;
    }

    template <class ScalarT, typename IdxT>
    ScalarT Bus<ScalarT, IdxT>::evaluateObjective(
        [[maybe_unused]] const ScalarT* local_variables) const
    {
      return ScalarT{0};
    }

    template <class ScalarT, typename IdxT>
    void Bus<ScalarT, IdxT>::evaluateConstraints(
        [[maybe_unused]] const ScalarT* local_variables,
        [[maybe_unused]] ScalarT*       local_values) const
    {
    }
  } // namespace OPF
} // namespace GridKit
