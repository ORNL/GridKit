#include "Shunt.hpp"

#include <cmath>
#include <string>

namespace GridKit
{
  namespace OPF
  {
    template <class ScalarT, typename IdxT>
    Shunt<ScalarT, IdxT>::Shunt(const DataT& data, const StateT& state)
      : data_(data),
        state_(state)
    {
    }

    template <class ScalarT, typename IdxT>
    typename Shunt<ScalarT, IdxT>::VariablesT& Shunt<ScalarT, IdxT>::variables()
    {
      return variables_;
    }

    template <class ScalarT, typename IdxT>
    const typename Shunt<ScalarT, IdxT>::VariablesT& Shunt<ScalarT, IdxT>::variables() const
    {
      return variables_;
    }

    template <class ScalarT, typename IdxT>
    typename Shunt<ScalarT, IdxT>::ConstraintsT& Shunt<ScalarT, IdxT>::constraints()
    {
      return constraints_;
    }

    template <class ScalarT, typename IdxT>
    const typename Shunt<ScalarT, IdxT>::ConstraintsT& Shunt<ScalarT, IdxT>::constraints() const
    {
      return constraints_;
    }

    template <class ScalarT, typename IdxT>
    IdxT Shunt<ScalarT, IdxT>::sizeInternalVariables() const
    {
      return static_cast<IdxT>(VariablesT::sizeInternal());
    }

    template <class ScalarT, typename IdxT>
    IdxT Shunt<ScalarT, IdxT>::sizeInternalConstraints() const
    {
      return static_cast<IdxT>(ConstraintsT::sizeInternal());
    }

    template <class ScalarT, typename IdxT>
    IdxT Shunt<ScalarT, IdxT>::nnz() const
    {
      return IdxT{2};
    }

    template <class ScalarT, typename IdxT>
    void Shunt<ScalarT, IdxT>::addJacobianSparsity(
        std::vector<typename ComponentT::JacobianEntry>& entries) const
    {
      const IdxT vm = variables_.template externalIndex<ShuntExternalVariables::VM>();
      entries.emplace_back(
          constraints_.template externalIndex<ShuntExternalConstraints::DIVP>(), vm);
      entries.emplace_back(
          constraints_.template externalIndex<ShuntExternalConstraints::DIVQ>(), vm);
    }

    template <class ScalarT, typename IdxT>
    int Shunt<ScalarT, IdxT>::setJacobianSlots(const std::vector<IdxT>& slots)
    {
      if (slots.size() != static_cast<std::size_t>(nnz()))
      {
        return 1;
      }
      jacobian_slots_ = slots;
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Shunt<ScalarT, IdxT>::addConstraints(const ScalarT* values,
                                             ScalarT*       constraints) const
    {
      if (!state_.online.value_or(true))
      {
        return 0;
      }

      const ScalarT vm  = variables_.template external<ShuntExternalVariables::VM>(values);
      const ScalarT vm2 = vm * vm;
      const ScalarT p   = -static_cast<ScalarT>(data_.G) * vm2;
      const ScalarT q   = static_cast<ScalarT>(data_.B) * vm2;
      if (!std::isfinite(static_cast<RealT>(p))
          || !std::isfinite(static_cast<RealT>(q)))
      {
        return 1;
      }

      constraints[constraints_.template externalIndex<ShuntExternalConstraints::DIVP>()] += p;
      constraints[constraints_.template externalIndex<ShuntExternalConstraints::DIVQ>()] += q;
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Shunt<ScalarT, IdxT>::addJacobian(const ScalarT* values,
                                          RealT*         jacobian_values) const
    {
      if (jacobian_slots_.size() != static_cast<std::size_t>(nnz()))
      {
        return 1;
      }
      if (!state_.online.value_or(true))
      {
        return 0;
      }

      const ScalarT vm           = variables_.template external<ShuntExternalVariables::VM>(values);
      const ScalarT derivative_p = -ScalarT{2}
                                   * static_cast<ScalarT>(data_.G) * vm;
      const ScalarT derivative_q = ScalarT{2}
                                   * static_cast<ScalarT>(data_.B) * vm;
      if (!std::isfinite(static_cast<RealT>(derivative_p))
          || !std::isfinite(static_cast<RealT>(derivative_q)))
      {
        return 1;
      }

      jacobian_values[jacobian_slots_[0]] += static_cast<RealT>(derivative_p);
      jacobian_values[jacobian_slots_[1]] += static_cast<RealT>(derivative_q);
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Shunt<ScalarT, IdxT>::updateSolutionState(const ScalarT*    values,
                                                  Model::StateData& state) const
    {
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

      const ScalarT vm = variables_.template external<ShuntExternalVariables::VM>(values);
      const ScalarT va = variables_.template external<ShuntExternalVariables::VA>(values);
      if (!std::isfinite(static_cast<RealT>(vm))
          || !std::isfinite(static_cast<RealT>(va)))
      {
        return 1;
      }

      const ScalarT vr        = vm * std::cos(va);
      const ScalarT vi        = vm * std::sin(va);
      const ScalarT current_r = -(static_cast<ScalarT>(data_.G) * vr
                                  - static_cast<ScalarT>(data_.B) * vi);
      const ScalarT current_i = -(static_cast<ScalarT>(data_.B) * vr
                                  + static_cast<ScalarT>(data_.G) * vi);
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

    template class Shunt<double, int>;
    template class Shunt<double, long int>;
    template class Shunt<double, std::size_t>;

  } // namespace OPF
} // namespace GridKit
