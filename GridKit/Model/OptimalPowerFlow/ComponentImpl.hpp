/**
 * @file ComponentImpl.hpp
 * @brief Definition of `Component`, included by model translation units.
 */

#pragma once

#include <algorithm>
#include <new>
#include <utility>
#include <vector>

#include <GridKit/Model/OptimalPowerFlow/Component.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    template <typename scalar_type, typename index_type>
    template <ComponentKernels<scalar_type, index_type> ModelT>
    Component<scalar_type, index_type>::Component(std::type_identity<ModelT>, std::string id)
      : id_(std::move(id)),
        size_{ModelT::NUM_INTERNAL_VARIABLES, ModelT::NUM_INTERNAL_EQUALITIES, ModelT::NUM_INTERNAL_INEQUALITIES}
    {
      for (IdxT port = 0; port < ModelT::NUM_PORTS; ++port)
      {
        ports_.push_back({{ModelT::NUM_PORT_VARIABLES[port], ModelT::NUM_PORT_EQUALITIES[port]}});
      }

      constexpr Block local = localSize<ModelT>();
      if (local_x_.resize(local.variables) != 0
          || local_c_.resize(local.equalities) != 0
          || equality_entries_.resize(local.equalities * local.variables) != 0
          || inequality_entries_.resize(local.inequalities * local.variables) != 0
          || hessian_entries_.resize(local.variables * local.variables) != 0)
      {
        throw std::bad_alloc{};
      }
    }

    /**
     * @brief Read the variables and contribute to the equality rows of
     * `component` from `offset` in its internal blocks
     */
    template <typename scalar_type, typename index_type>
    int Component<scalar_type, index_type>::connect(IdxT port, Component& component, const Block& offset)
    {
      if (bound_ || port >= ports_.size())
      {
        return 1;
      }

      const Block& size = ports_[port].size;
      const auto   fits = [](IdxT start, IdxT count, IdxT available)
      { return start <= available && count <= available - start; };
      if (!fits(offset.variables, size.variables, component.size_.variables)
          || !fits(offset.equalities, size.equalities, component.size_.equalities)
          || offset.inequalities != 0)
      {
        return 1;
      }

      ports_[port].component = &component;
      ports_[port].offset    = offset;
      return 0;
    }

    /// Local variables: internal, then those each port reads
    template <typename scalar_type, typename index_type>
    void Component<scalar_type, index_type>::gather(ScalarT* x) const
    {
      x = std::copy_n(this->x_.getData(), size_.variables, x);
      for (const Port& port : ports_)
      {
        x = std::copy_n(port.component->x_.getData() + port.offset.variables, port.size.variables, x);
      }
    }

    /**
     * @brief Global index of each local entry of `block`: internal, then those
     * of each port
     */
    template <typename scalar_type, typename index_type>
    void Component<scalar_type, index_type>::globalIndices(IdxT Block::* block, IdxT* indices) const
    {
      IdxT local = 0;
      for (IdxT i = 0; i < size_.*block; ++i)
      {
        indices[local++] = offset_.*block + i;
      }
      for (const Port& port : ports_)
      {
        for (IdxT i = 0; i < port.size.*block; ++i)
        {
          indices[local++] = port.component->offset_.*block + port.offset.*block + i;
        }
      }
    }

    /**
     * @brief Add the local entries of `block` to `vector`: internal, then
     * those of each port
     */
    template <typename scalar_type, typename index_type>
    template <typename ValueT, typename MemberT>
    void Component<scalar_type, index_type>::scatter(IdxT Block::* block, MemberT Component::* vector, const ValueT* local)
    {
      ValueT* values = (this->*vector).getData();
      for (IdxT i = 0; i < size_.*block; ++i)
      {
        values[i] += *local++;
      }
      for (const Port& port : ports_)
      {
        ValueT* target = (port.component->*vector).getData() + port.offset.*block;
        for (IdxT i = 0; i < port.size.*block; ++i)
        {
          target[i] += *local++;
        }
      }
    }

    template <typename scalar_type, typename index_type>
    int Component<scalar_type, index_type>::evaluateObjective()
    {
      gather(local_x_.getData());
      this->f_ = objective(local_x_.getData());
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Component<scalar_type, index_type>::evaluateEqualities()
    {
      gather(local_x_.getData());
      local_c_.setToZero();
      internalEqualities(local_x_.getData(), local_c_.getData());
      externalEqualities(local_x_.getData(), local_c_.getData() + size_.equalities);
      scatter(&Block::equalities, &Component::c_, local_c_.getData());
      return 0;
    }

    /**
     * @brief Inequality rows are all internal, so the kernel writes the view
     */
    template <typename scalar_type, typename index_type>
    int Component<scalar_type, index_type>::evaluateInequalities()
    {
      gather(local_x_.getData());
      internalInequalities(local_x_.getData(), this->d_.getData());
      return 0;
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
