#include "Generator.hpp"

#include <cmath>
#include <limits>
#include <string>

namespace GridKit
{
  namespace OPF
  {
    template <class ScalarT, typename IdxT>
    Generator<ScalarT, IdxT>::Generator(const DataT& data, const StateT& state)
      : data_(data),
        state_(state)
    {
    }

    template <class ScalarT, typename IdxT>
    typename Generator<ScalarT, IdxT>::VariablesT& Generator<ScalarT, IdxT>::variables()
    {
      return variables_;
    }

    template <class ScalarT, typename IdxT>
    const typename Generator<ScalarT, IdxT>::VariablesT& Generator<ScalarT, IdxT>::variables() const
    {
      return variables_;
    }

    template <class ScalarT, typename IdxT>
    typename Generator<ScalarT, IdxT>::ConstraintsT& Generator<ScalarT, IdxT>::constraints()
    {
      return constraints_;
    }

    template <class ScalarT, typename IdxT>
    const typename Generator<ScalarT, IdxT>::ConstraintsT& Generator<ScalarT, IdxT>::constraints() const
    {
      return constraints_;
    }

    template <class ScalarT, typename IdxT>
    IdxT Generator<ScalarT, IdxT>::sizeInternalVariables() const
    {
      return static_cast<IdxT>(VariablesT::sizeInternal());
    }

    template <class ScalarT, typename IdxT>
    IdxT Generator<ScalarT, IdxT>::sizeInternalConstraints() const
    {
      return static_cast<IdxT>(ConstraintsT::sizeInternal());
    }

    template <class ScalarT, typename IdxT>
    IdxT Generator<ScalarT, IdxT>::nnz() const
    {
      return IdxT{2};
    }

    template <class ScalarT, typename IdxT>
    void Generator<ScalarT, IdxT>::setVariableOffset(IdxT offset)
    {
      variables_.setInternalOffset(offset);
    }

    template <class ScalarT, typename IdxT>
    void Generator<ScalarT, IdxT>::addJacobianSparsity(
        std::vector<typename ComponentT::JacobianEntry>& entries) const
    {
      entries.emplace_back(
          constraints_.template externalIndex<GeneratorExternalConstraints::DIVP>(),
          variables_.template internalIndex<GeneratorInternalVariables::PG>());
      entries.emplace_back(
          constraints_.template externalIndex<GeneratorExternalConstraints::DIVQ>(),
          variables_.template internalIndex<GeneratorInternalVariables::QG>());
    }

    template <class ScalarT, typename IdxT>
    int Generator<ScalarT, IdxT>::setJacobianSlots(const std::vector<IdxT>& slots)
    {
      if (slots.size() != static_cast<std::size_t>(nnz()))
      {
        return 1;
      }
      jacobian_slots_ = slots;
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Generator<ScalarT, IdxT>::initialize(ScalarT* values) const
    {
      if ((state_.p.has_value() && !std::isfinite(*state_.p))
          || (state_.q.has_value() && !std::isfinite(*state_.q)))
      {
        return 1;
      }

      const bool online                                                    = state_.online.value_or(true);
      variables_.template internal<GeneratorInternalVariables::PG>(values) = online ? static_cast<ScalarT>(state_.p.value_or(0.0)) : ScalarT{0};
      variables_.template internal<GeneratorInternalVariables::QG>(values) = online ? static_cast<ScalarT>(state_.q.value_or(0.0)) : ScalarT{0};
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Generator<ScalarT, IdxT>::setVariableBounds(RealT* lower,
                                                    RealT* upper) const
    {
      const IdxT pg = variables_.template internalIndex<GeneratorInternalVariables::PG>();
      const IdxT qg = variables_.template internalIndex<GeneratorInternalVariables::QG>();
      if (!state_.online.value_or(true))
      {
        lower[pg] = RealT{0};
        upper[pg] = RealT{0};
        lower[qg] = RealT{0};
        upper[qg] = RealT{0};
        return 0;
      }

      const RealT infinity = std::numeric_limits<RealT>::infinity();
      const RealT pg_lower = data_.pmin.value_or(-infinity);
      const RealT pg_upper = data_.pmax.value_or(infinity);
      const RealT qg_lower = data_.qmin.value_or(-infinity);
      const RealT qg_upper = data_.qmax.value_or(infinity);
      if ((data_.pmin.has_value() && !std::isfinite(*data_.pmin))
          || (data_.pmax.has_value() && !std::isfinite(*data_.pmax))
          || (data_.qmin.has_value() && !std::isfinite(*data_.qmin))
          || (data_.qmax.has_value() && !std::isfinite(*data_.qmax))
          || pg_lower > pg_upper || qg_lower > qg_upper)
      {
        return 1;
      }

      lower[pg] = pg_lower;
      upper[pg] = pg_upper;
      lower[qg] = qg_lower;
      upper[qg] = qg_upper;
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Generator<ScalarT, IdxT>::addObjective(const ScalarT* values,
                                               RealT&         objective) const
    {
      if (!state_.online.value_or(true))
      {
        return 0;
      }
      if (!std::isfinite(data_.c0) || !std::isfinite(data_.c1)
          || !std::isfinite(data_.c2))
      {
        return 1;
      }

      const ScalarT pg           = variables_.template internal<GeneratorInternalVariables::PG>(values);
      const ScalarT contribution = static_cast<ScalarT>(data_.c0)
                                   + static_cast<ScalarT>(data_.c1) * pg
                                   + static_cast<ScalarT>(data_.c2) * pg * pg;
      const RealT updated_objective = objective + static_cast<RealT>(contribution);
      if (!std::isfinite(static_cast<RealT>(pg))
          || !std::isfinite(static_cast<RealT>(contribution))
          || !std::isfinite(updated_objective))
      {
        return 1;
      }
      objective = updated_objective;
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Generator<ScalarT, IdxT>::addObjectiveGradient(
        const ScalarT* values,
        ScalarT*       gradient) const
    {
      if (!state_.online.value_or(true))
      {
        return 0;
      }
      if (!std::isfinite(data_.c1) || !std::isfinite(data_.c2))
      {
        return 1;
      }

      const IdxT    pg_index     = variables_.template internalIndex<GeneratorInternalVariables::PG>();
      const ScalarT pg           = values[pg_index];
      const ScalarT contribution = static_cast<ScalarT>(data_.c1)
                                   + ScalarT{2}
                                         * static_cast<ScalarT>(data_.c2) * pg;
      const ScalarT updated_gradient = gradient[pg_index] + contribution;
      if (!std::isfinite(static_cast<RealT>(pg))
          || !std::isfinite(static_cast<RealT>(contribution))
          || !std::isfinite(static_cast<RealT>(updated_gradient)))
      {
        return 1;
      }
      gradient[pg_index] = updated_gradient;
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Generator<ScalarT, IdxT>::addConstraints(const ScalarT* values,
                                                 ScalarT*       constraints) const
    {
      constraints[constraints_.template externalIndex<GeneratorExternalConstraints::DIVP>()] += variables_.template internal<GeneratorInternalVariables::PG>(values);
      constraints[constraints_.template externalIndex<GeneratorExternalConstraints::DIVQ>()] += variables_.template internal<GeneratorInternalVariables::QG>(values);
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Generator<ScalarT, IdxT>::addJacobian(
        [[maybe_unused]] const ScalarT* values,
        RealT*                          jacobian_values) const
    {
      if (jacobian_slots_.size() != static_cast<std::size_t>(nnz()))
      {
        return 1;
      }

      jacobian_values[jacobian_slots_[0]] += RealT{1};
      jacobian_values[jacobian_slots_[1]] += RealT{1};
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Generator<ScalarT, IdxT>::updateSolutionState(
        const ScalarT*    values,
        Model::StateData& state) const
    {
      const ScalarT pg = variables_.template internal<GeneratorInternalVariables::PG>(values);
      const ScalarT qg = variables_.template internal<GeneratorInternalVariables::QG>(values);
      if (!std::isfinite(static_cast<RealT>(pg))
          || !std::isfinite(static_cast<RealT>(qg)))
      {
        return 1;
      }

      const auto bus_entry = state.buses.find("bus_id_" + std::to_string(data_.bus));
      if (bus_entry == state.buses.end())
      {
        return 1;
      }

      if (!state_.online.value_or(true))
      {
        auto& device = state.devices[data_.id];
        device.p     = static_cast<double>(pg);
        device.q     = static_cast<double>(qg);

        auto& injection = bus_entry->second.injections[data_.id];
        injection.ir    = 0.0;
        injection.ii    = 0.0;
        return 0;
      }

      const ScalarT vm = variables_.template external<GeneratorExternalVariables::VM>(values);
      const ScalarT va = variables_.template external<GeneratorExternalVariables::VA>(values);
      if (!(vm > ScalarT{0}) || !std::isfinite(static_cast<RealT>(vm))
          || !std::isfinite(static_cast<RealT>(va)))
      {
        return 1;
      }

      const ScalarT vr        = vm * std::cos(va);
      const ScalarT vi        = vm * std::sin(va);
      const ScalarT vm2       = vm * vm;
      const ScalarT current_r = (pg * vr + qg * vi) / vm2;
      const ScalarT current_i = (pg * vi - qg * vr) / vm2;
      if (!std::isfinite(static_cast<RealT>(current_r))
          || !std::isfinite(static_cast<RealT>(current_i)))
      {
        return 1;
      }

      auto& device = state.devices[data_.id];
      device.p     = static_cast<double>(pg);
      device.q     = static_cast<double>(qg);

      auto& injection = bus_entry->second.injections[data_.id];
      injection.ir    = static_cast<double>(current_r);
      injection.ii    = static_cast<double>(current_i);
      return 0;
    }

    template class Generator<double, int>;
    template class Generator<double, long int>;
    template class Generator<double, std::size_t>;

  } // namespace OPF
} // namespace GridKit
