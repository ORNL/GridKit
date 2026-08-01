#include "Load.hpp"

#include <cmath>
#include <string>

namespace GridKit
{
  namespace OPF
  {
    template <class ScalarT, typename IdxT>
    Load<ScalarT, IdxT>::Load(const DataT& data, const StateT& state)
      : data_(data),
        state_(state)
    {
    }

    template <class ScalarT, typename IdxT>
    typename Load<ScalarT, IdxT>::VariablesT& Load<ScalarT, IdxT>::variables()
    {
      return variables_;
    }

    template <class ScalarT, typename IdxT>
    const typename Load<ScalarT, IdxT>::VariablesT& Load<ScalarT, IdxT>::variables() const
    {
      return variables_;
    }

    template <class ScalarT, typename IdxT>
    typename Load<ScalarT, IdxT>::ConstraintsT& Load<ScalarT, IdxT>::constraints()
    {
      return constraints_;
    }

    template <class ScalarT, typename IdxT>
    const typename Load<ScalarT, IdxT>::ConstraintsT& Load<ScalarT, IdxT>::constraints() const
    {
      return constraints_;
    }

    template <class ScalarT, typename IdxT>
    IdxT Load<ScalarT, IdxT>::sizeInternalVariables() const
    {
      return static_cast<IdxT>(VariablesT::sizeInternal());
    }

    template <class ScalarT, typename IdxT>
    IdxT Load<ScalarT, IdxT>::sizeInternalConstraints() const
    {
      return static_cast<IdxT>(ConstraintsT::sizeInternal());
    }

    template <class ScalarT, typename IdxT>
    IdxT Load<ScalarT, IdxT>::nnz() const
    {
      return IdxT{0};
    }

    template <class ScalarT, typename IdxT>
    int Load<ScalarT, IdxT>::addConstraints(
        [[maybe_unused]] const ScalarT* values,
        ScalarT*                        constraints) const
    {
      if (!state_.p.has_value() || !state_.q.has_value()
          || !std::isfinite(*state_.p) || !std::isfinite(*state_.q))
      {
        return 1;
      }

      const ScalarT online                                                               = state_.online.value_or(true) ? ScalarT{1} : ScalarT{0};
      constraints[constraints_.template externalIndex<LoadExternalConstraints::DIVP>()] += online * static_cast<ScalarT>(*state_.p);
      constraints[constraints_.template externalIndex<LoadExternalConstraints::DIVQ>()] += online * static_cast<ScalarT>(*state_.q);
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Load<ScalarT, IdxT>::updateSolutionState(const ScalarT*    values,
                                                 Model::StateData& state) const
    {
      if (!state_.p.has_value() || !state_.q.has_value()
          || !std::isfinite(*state_.p) || !std::isfinite(*state_.q))
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
        auto& injection = bus_entry->second.injections[data_.id];
        injection.ir    = 0.0;
        injection.ii    = 0.0;
        return 0;
      }

      const ScalarT vm = variables_.template external<LoadExternalVariables::VM>(values);
      const ScalarT va = variables_.template external<LoadExternalVariables::VA>(values);
      if (!(vm > ScalarT{0}) || !std::isfinite(static_cast<RealT>(vm))
          || !std::isfinite(static_cast<RealT>(va)))
      {
        return 1;
      }

      const ScalarT p         = static_cast<ScalarT>(*state_.p);
      const ScalarT q         = static_cast<ScalarT>(*state_.q);
      const ScalarT cos_va    = std::cos(va);
      const ScalarT sin_va    = std::sin(va);
      const ScalarT current_r = (p * cos_va + q * sin_va) / vm;
      const ScalarT current_i = (p * sin_va - q * cos_va) / vm;
      if (!std::isfinite(static_cast<RealT>(current_r))
          || !std::isfinite(static_cast<RealT>(current_i)))
      {
        return 1;
      }

      auto& injection = bus_entry->second.injections[data_.id];
      injection.ir    = static_cast<double>(current_r);
      injection.ii    = static_cast<double>(current_i);
      return 0;
    }

    template class Load<double, int>;
    template class Load<double, long int>;
    template class Load<double, std::size_t>;

  } // namespace OPF
} // namespace GridKit
