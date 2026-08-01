#include "Bus.hpp"

#include <cmath>
#include <limits>
#include <string>

namespace GridKit
{
  namespace OPF
  {
    template <class ScalarT, typename IdxT>
    Bus<ScalarT, IdxT>::Bus(const DataT& data, const StateT& state)
      : data_(data),
        state_(state)
    {
    }

    template <class ScalarT, typename IdxT>
    typename Bus<ScalarT, IdxT>::VariablesT& Bus<ScalarT, IdxT>::variables()
    {
      return variables_;
    }

    template <class ScalarT, typename IdxT>
    const typename Bus<ScalarT, IdxT>::VariablesT& Bus<ScalarT, IdxT>::variables() const
    {
      return variables_;
    }

    template <class ScalarT, typename IdxT>
    typename Bus<ScalarT, IdxT>::ConstraintsT& Bus<ScalarT, IdxT>::constraints()
    {
      return constraints_;
    }

    template <class ScalarT, typename IdxT>
    const typename Bus<ScalarT, IdxT>::ConstraintsT& Bus<ScalarT, IdxT>::constraints() const
    {
      return constraints_;
    }

    template <class ScalarT, typename IdxT>
    IdxT Bus<ScalarT, IdxT>::sizeInternalVariables() const
    {
      return static_cast<IdxT>(VariablesT::sizeInternal());
    }

    template <class ScalarT, typename IdxT>
    IdxT Bus<ScalarT, IdxT>::sizeInternalConstraints() const
    {
      return static_cast<IdxT>(ConstraintsT::sizeInternal());
    }

    template <class ScalarT, typename IdxT>
    IdxT Bus<ScalarT, IdxT>::nnz() const
    {
      return IdxT{0};
    }

    template <class ScalarT, typename IdxT>
    void Bus<ScalarT, IdxT>::setVariableOffset(IdxT offset)
    {
      variables_.setInternalOffset(offset);
    }

    template <class ScalarT, typename IdxT>
    void Bus<ScalarT, IdxT>::setConstraintOffset(IdxT offset)
    {
      constraints_.setInternalOffset(offset);
    }

    template <class ScalarT, typename IdxT>
    int Bus<ScalarT, IdxT>::initialize(ScalarT* values) const
    {
      if (!state_.vr.has_value() || !state_.vi.has_value()
          || !std::isfinite(*state_.vr) || !std::isfinite(*state_.vi))
      {
        return 1;
      }

      const RealT vm = std::hypot(static_cast<RealT>(*state_.vr),
                                  static_cast<RealT>(*state_.vi));
      if (!(vm > RealT{0}) || !std::isfinite(vm))
      {
        return 1;
      }

      variables_.template internal<BusInternalVariables::VM>(values) = static_cast<ScalarT>(vm);
      variables_.template internal<BusInternalVariables::VA>(values) = static_cast<ScalarT>(std::atan2(*state_.vi, *state_.vr));
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Bus<ScalarT, IdxT>::setVariableBounds(RealT* lower,
                                              RealT* upper) const
    {
      const RealT infinity = std::numeric_limits<RealT>::infinity();
      const RealT vm_lower = data_.vmin.value_or(-infinity);
      const RealT vm_upper = data_.vmax.value_or(infinity);
      if (!std::isfinite(data_.kv) || !(data_.kv > RealT{0})
          || (data_.vmin.has_value()
              && (!std::isfinite(*data_.vmin) || *data_.vmin < RealT{0}))
          || (data_.vmax.has_value()
              && (!std::isfinite(*data_.vmax) || !(*data_.vmax > RealT{0})))
          || vm_lower > vm_upper)
      {
        return 1;
      }

      const IdxT vm = variables_.template internalIndex<BusInternalVariables::VM>();
      const IdxT va = variables_.template internalIndex<BusInternalVariables::VA>();
      lower[vm]     = vm_lower;
      upper[vm]     = vm_upper;
      lower[va]     = -infinity;
      upper[va]     = infinity;
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Bus<ScalarT, IdxT>::updateSolutionState(const ScalarT*    values,
                                                Model::StateData& state) const
    {
      const ScalarT vm = variables_.template internal<BusInternalVariables::VM>(values);
      const ScalarT va = variables_.template internal<BusInternalVariables::VA>(values);
      if (!std::isfinite(static_cast<RealT>(vm))
          || !std::isfinite(static_cast<RealT>(va)))
      {
        return 1;
      }

      const auto bus_entry = state.buses.find("bus_id_" + std::to_string(data_.number));
      if (bus_entry == state.buses.end())
      {
        return 1;
      }

      bus_entry->second.vr = static_cast<double>(vm * std::cos(va));
      bus_entry->second.vi = static_cast<double>(vm * std::sin(va));
      return 0;
    }

    template class Bus<double, int>;
    template class Bus<double, long int>;
    template class Bus<double, std::size_t>;

  } // namespace OPF
} // namespace GridKit
