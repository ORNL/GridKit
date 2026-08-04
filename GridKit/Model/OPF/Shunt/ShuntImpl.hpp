#pragma once

#include <cmath>
#include <string>

#include "Shunt.hpp"

namespace GridKit
{
  namespace OPF
  {
    template <class ScalarT, typename IdxT>
    Shunt<ScalarT, IdxT>::Shunt(
        const DataT&                       data,
        const StateT&                      state,
        IdxT                               voltage_magnitude_index,
        std::array<IdxT, CONSTRAINT_COUNT> constraint_indices,
        std::array<IdxT, 2>                bus_voltage_indices)
      : data_(data),
        in_service_(state.online.value_or(true) ? RealT{1} : RealT{0}),
        variable_indices_{voltage_magnitude_index},
        constraint_indices_(constraint_indices),
        bus_voltage_indices_(bus_voltage_indices)
    {
    }

    template <class ScalarT, typename IdxT>
    std::span<const IdxT> Shunt<ScalarT, IdxT>::variableIndices() const
    {
      return variable_indices_;
    }

    template <class ScalarT, typename IdxT>
    std::span<const IdxT> Shunt<ScalarT, IdxT>::constraintIndices() const
    {
      return constraint_indices_;
    }

    template <class ScalarT, typename IdxT>
    std::span<const typename Shunt<ScalarT, IdxT>::JacobianEntryT>
    Shunt<ScalarT, IdxT>::jacobianEntries() const
    {
      return JACOBIAN_ENTRIES;
    }

    template <class ScalarT, typename IdxT>
    std::span<const typename Shunt<ScalarT, IdxT>::HessianEntryT>
    Shunt<ScalarT, IdxT>::hessianEntries() const
    {
      return HESSIAN_ENTRIES;
    }

    template <class ScalarT, typename IdxT>
    int Shunt<ScalarT, IdxT>::gatherVariables(const ScalarT* global_variables)
    {
      if (global_variables == nullptr)
      {
        return 1;
      }
      variables_[0] = global_variables[variable_indices_[0]];
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Shunt<ScalarT, IdxT>::evaluateObjective()
    {
      objective_ = RealT{0};
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Shunt<ScalarT, IdxT>::evaluateObjectiveGradient()
    {
      objective_gradient_values_.fill(ScalarT{0});
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Shunt<ScalarT, IdxT>::evaluateConstraints()
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
    typename Shunt<ScalarT, IdxT>::RealT Shunt<ScalarT, IdxT>::objective() const
    {
      return objective_;
    }

    template <class ScalarT, typename IdxT>
    std::span<const ScalarT> Shunt<ScalarT, IdxT>::objectiveGradientValues() const
    {
      return objective_gradient_values_;
    }

    template <class ScalarT, typename IdxT>
    std::span<const ScalarT> Shunt<ScalarT, IdxT>::constraintValues() const
    {
      return constraint_values_;
    }

    template <class ScalarT, typename IdxT>
    std::span<const typename Shunt<ScalarT, IdxT>::RealT>
    Shunt<ScalarT, IdxT>::jacobianValues() const
    {
      return jacobian_values_;
    }

    template <class ScalarT, typename IdxT>
    std::span<const typename Shunt<ScalarT, IdxT>::RealT>
    Shunt<ScalarT, IdxT>::hessianValues() const
    {
      return hessian_values_;
    }

    template <class ScalarT, typename IdxT>
    bool Shunt<ScalarT, IdxT>::hasJacobian() const
    {
      return true;
    }

    template <class ScalarT, typename IdxT>
    bool Shunt<ScalarT, IdxT>::hasHessian() const
    {
      return true;
    }

    template <class ScalarT, typename IdxT>
    int Shunt<ScalarT, IdxT>::updateSolutionState(
        const ScalarT*    global_variables,
        Model::StateData& state) const
    {
      if (global_variables == nullptr)
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
        const ScalarT vr        = magnitude * std::cos(angle);
        const ScalarT vi        = magnitude * std::sin(angle);
        current_r               = -(static_cast<ScalarT>(data_.G) * vr
                      - static_cast<ScalarT>(data_.B) * vi);
        current_i               = -(static_cast<ScalarT>(data_.B) * vr
                      + static_cast<ScalarT>(data_.G) * vi);
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
    ScalarT Shunt<ScalarT, IdxT>::evaluateObjective(
        [[maybe_unused]] const ScalarT* local_variables) const
    {
      return ScalarT{0};
    }

    template <class ScalarT, typename IdxT>
    void Shunt<ScalarT, IdxT>::evaluateConstraints(
        const ScalarT* local_variables,
        ScalarT*       local_values) const
    {
      const ScalarT v2 = local_variables[0] * local_variables[0];
      local_values[0]  = -static_cast<ScalarT>(in_service_ * data_.G) * v2;
      local_values[1]  = static_cast<ScalarT>(in_service_ * data_.B) * v2;
    }
  } // namespace OPF
} // namespace GridKit
