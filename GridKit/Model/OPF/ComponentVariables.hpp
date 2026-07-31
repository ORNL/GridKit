#pragma once

#include <array>
#include <cstddef>
#include <stdexcept>
#include <type_traits>

#include <GridKit/Constants.hpp>

namespace GridKit
{
  namespace OPF
  {
    enum class NoVariables : std::size_t
    {
      MAXIMUM
    };

    template <typename EnumT>
    concept VariableEnum = std::is_enum_v<EnumT> && std::is_same_v<std::underlying_type_t<EnumT>, std::size_t> && requires { EnumT::MAXIMUM; };

    /**
     * @brief Global variable indices owned and consumed by an OPF component.
     *
     * Internal variables are entries in the system decision vector owned by
     * the component. External variables are entries owned by another component
     * and connected by System while composing the network.
     */
    template <class scalar_type,
              typename index_type,
              VariableEnum InternalVariables,
              VariableEnum ExternalVariables>
    class ComponentVariables
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;

      static constexpr std::size_t sizeInternal()
      {
        return static_cast<std::size_t>(InternalVariables::MAXIMUM);
      }

      static constexpr std::size_t sizeExternal()
      {
        return static_cast<std::size_t>(ExternalVariables::MAXIMUM);
      }

      ComponentVariables()
      {
        internal_indices_.fill(GridKit::INVALID_INDEX<IdxT>);
        external_indices_.fill(GridKit::INVALID_INDEX<IdxT>);
      }

      void setInternalOffset(IdxT offset)
      {
        for (std::size_t i = 0; i < sizeInternal(); ++i)
        {
          internal_indices_[i] = offset + static_cast<IdxT>(i);
        }
      }

      template <ExternalVariables variable>
      void bindExternal(IdxT index)
      {
        static_assert(variable < ExternalVariables::MAXIMUM);
        external_indices_[static_cast<std::size_t>(variable)] = index;
      }

      template <InternalVariables variable>
      IdxT internalIndex() const
      {
        static_assert(variable < InternalVariables::MAXIMUM);
        return checkedIndex(internal_indices_[static_cast<std::size_t>(variable)]);
      }

      template <ExternalVariables variable>
      IdxT externalIndex() const
      {
        static_assert(variable < ExternalVariables::MAXIMUM);
        return checkedIndex(external_indices_[static_cast<std::size_t>(variable)]);
      }

      template <InternalVariables variable>
      ScalarT& internal(ScalarT* values) const
      {
        return values[internalIndex<variable>()];
      }

      template <InternalVariables variable>
      const ScalarT& internal(const ScalarT* values) const
      {
        return values[internalIndex<variable>()];
      }

      template <ExternalVariables variable>
      const ScalarT& external(const ScalarT* values) const
      {
        return values[externalIndex<variable>()];
      }

    private:
      static IdxT checkedIndex(IdxT index)
      {
        if (index == GridKit::INVALID_INDEX<IdxT>)
        {
          throw std::logic_error("OPF component variable is not bound");
        }
        return index;
      }

      std::array<IdxT, sizeInternal()> internal_indices_;
      std::array<IdxT, sizeExternal()> external_indices_;
    };

  } // namespace OPF
} // namespace GridKit
