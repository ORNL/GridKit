#pragma once

#include <cmath>
#include <limits>
#include <string>

#include "Load.hpp"

namespace GridKit
{
  namespace OPF
  {
    template <class ScalarT, typename IdxT>
    Load<ScalarT, IdxT>::Load(
        const DataT&                       data,
        const StateT&                      state,
        std::array<IdxT, CONSTRAINT_COUNT> constraint_indices,
        std::array<IdxT, 2>                bus_voltage_indices)
      : data_(data),
        in_service_(state.online.value_or(true) ? RealT{1} : RealT{0}),
        p_(state.p.value_or(std::numeric_limits<RealT>::quiet_NaN())),
        q_(state.q.value_or(std::numeric_limits<RealT>::quiet_NaN())),
        constraint_indices_(constraint_indices),
        bus_voltage_indices_(bus_voltage_indices)
    {
    }

    template <class ScalarT, typename IdxT>
    std::span<const IdxT> Load<ScalarT, IdxT>::variableIndices() const
    {
      return {};
    }

    template <class ScalarT, typename IdxT>
    std::span<const IdxT> Load<ScalarT, IdxT>::constraintIndices() const
    {
      return constraint_indices_;
    }

    template <class ScalarT, typename IdxT>
    std::span<const typename Load<ScalarT, IdxT>::JacobianEntryT>
    Load<ScalarT, IdxT>::jacobianEntries() const
    {
      return {};
    }

    template <class ScalarT, typename IdxT>
    std::span<const typename Load<ScalarT, IdxT>::HessianEntryT>
    Load<ScalarT, IdxT>::hessianEntries() const
    {
      return {};
    }

    template <class ScalarT, typename IdxT>
    int Load<ScalarT, IdxT>::gatherVariables(
        [[maybe_unused]] const ScalarT* global_variables)
    {
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Load<ScalarT, IdxT>::evaluateObjective()
    {
      objective_ = RealT{0};
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Load<ScalarT, IdxT>::evaluateObjectiveGradient()
    {
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Load<ScalarT, IdxT>::evaluateConstraints()
    {
      evaluateConstraints(nullptr, constraint_values_.data());
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
    int Load<ScalarT, IdxT>::evaluateJacobian()
    {
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Load<ScalarT, IdxT>::evaluateHessian(
        [[maybe_unused]] RealT objective_factor,
        std::span<const RealT> local_multipliers)
    {
      return local_multipliers.size() == CONSTRAINT_COUNT ? 0 : 1;
    }

    template <class ScalarT, typename IdxT>
    typename Load<ScalarT, IdxT>::RealT Load<ScalarT, IdxT>::objective() const
    {
      return objective_;
    }

    template <class ScalarT, typename IdxT>
    std::span<const ScalarT> Load<ScalarT, IdxT>::objectiveGradientValues() const
    {
      return {};
    }

    template <class ScalarT, typename IdxT>
    std::span<const ScalarT> Load<ScalarT, IdxT>::constraintValues() const
    {
      return constraint_values_;
    }

    template <class ScalarT, typename IdxT>
    std::span<const typename Load<ScalarT, IdxT>::RealT>
    Load<ScalarT, IdxT>::jacobianValues() const
    {
      return {};
    }

    template <class ScalarT, typename IdxT>
    std::span<const typename Load<ScalarT, IdxT>::RealT>
    Load<ScalarT, IdxT>::hessianValues() const
    {
      return {};
    }

    template <class ScalarT, typename IdxT>
    bool Load<ScalarT, IdxT>::hasJacobian() const
    {
      return true;
    }

    template <class ScalarT, typename IdxT>
    bool Load<ScalarT, IdxT>::hasHessian() const
    {
      return true;
    }

    template <class ScalarT, typename IdxT>
    int Load<ScalarT, IdxT>::updateSolutionState(
        const ScalarT*    global_variables,
        Model::StateData& state) const
    {
      if (global_variables == nullptr || !std::isfinite(p_) || !std::isfinite(q_))
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
        current_r = (static_cast<ScalarT>(p_) * cosine
                     + static_cast<ScalarT>(q_) * sine)
                    / magnitude;
        current_i = (static_cast<ScalarT>(p_) * sine
                     - static_cast<ScalarT>(q_) * cosine)
                    / magnitude;
      }

      if (!std::isfinite(static_cast<RealT>(current_r))
          || !std::isfinite(static_cast<RealT>(current_i)))
      {
        return 1;
      }

      auto& injection = bus->second.injections[data_.id];
      injection.ir    = static_cast<double>(current_r);
      injection.ii    = static_cast<double>(current_i);
      return 0;
    }

    template <class ScalarT, typename IdxT>
    ScalarT Load<ScalarT, IdxT>::evaluateObjective(
        [[maybe_unused]] const ScalarT* local_variables) const
    {
      return ScalarT{0};
    }

    template <class ScalarT, typename IdxT>
    void Load<ScalarT, IdxT>::evaluateConstraints(
        [[maybe_unused]] const ScalarT* local_variables,
        ScalarT*                        local_values) const
    {
      local_values[0] = static_cast<ScalarT>(in_service_ * p_);
      local_values[1] = static_cast<ScalarT>(in_service_ * q_);
    }
  } // namespace OPF
} // namespace GridKit
