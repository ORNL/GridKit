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
    enum class NoConstraints : std::size_t
    {
      MAXIMUM
    };

    template <typename EnumT>
    concept ConstraintEnum = std::is_enum_v<EnumT> && std::is_same_v<std::underlying_type_t<EnumT>, std::size_t> && requires { EnumT::MAXIMUM; };

    /**
     * @brief Global constraint rows owned and consumed by an OPF component.
     *
     * Internal constraints are rows owned by the component. External
     * constraints are rows owned elsewhere, such as a bus power-balance row,
     * to which the component adds a contribution.
     */
    template <typename index_type,
              ConstraintEnum InternalConstraints,
              ConstraintEnum ExternalConstraints>
    class ComponentConstraints
    {
    public:
      using IdxT = index_type;

      static constexpr std::size_t sizeInternal()
      {
        return static_cast<std::size_t>(InternalConstraints::MAXIMUM);
      }

      static constexpr std::size_t sizeExternal()
      {
        return static_cast<std::size_t>(ExternalConstraints::MAXIMUM);
      }

      ComponentConstraints()
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

      template <ExternalConstraints constraint>
      void bindExternal(IdxT index)
      {
        static_assert(constraint < ExternalConstraints::MAXIMUM);
        external_indices_[static_cast<std::size_t>(constraint)] = index;
      }

      template <InternalConstraints constraint>
      IdxT internalIndex() const
      {
        static_assert(constraint < InternalConstraints::MAXIMUM);
        return checkedIndex(internal_indices_[static_cast<std::size_t>(constraint)]);
      }

      template <ExternalConstraints constraint>
      IdxT externalIndex() const
      {
        static_assert(constraint < ExternalConstraints::MAXIMUM);
        return checkedIndex(external_indices_[static_cast<std::size_t>(constraint)]);
      }

    private:
      static IdxT checkedIndex(IdxT index)
      {
        if (index == GridKit::INVALID_INDEX<IdxT>)
        {
          throw std::logic_error("OPF component constraint is not bound");
        }
        return index;
      }

      std::array<IdxT, sizeInternal()> internal_indices_;
      std::array<IdxT, sizeExternal()> external_indices_;
    };

  } // namespace OPF
} // namespace GridKit
