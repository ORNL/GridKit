/**
 * @file ComponentImpl.hpp
 * @brief Definition of `Component`, included by model translation units.
 */

#pragma once

#include <algorithm>
#include <new>
#include <numeric>
#include <utility>

#include <GridKit/Model/OptimalPowerFlow/Component.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    template <typename scalar_type, typename index_type>
    template <ComponentKernels<scalar_type, index_type> ModelT>
    Component<scalar_type, index_type>::Component(std::type_identity<ModelT>, std::string id)
      : id_(std::move(id)),
        size_(internalSize<ModelT>())
    {
      for (IdxT port = 0; port < ModelT::NUM_PORT_VARIABLES.size(); ++port)
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

    template <typename scalar_type, typename index_type>
    int Component<scalar_type, index_type>::evaluateObjective()
    {
      ScalarT* x = local_x_.getData();
      gather(x);
      f_ = objective(x, x + size_.variables);
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Component<scalar_type, index_type>::evaluateEqualities()
    {
      ScalarT* x = local_x_.getData();
      ScalarT* c = local_c_.getData();
      gather(x);
      local_c_.setToZero();
      internalEqualities(x, x + size_.variables, c);
      externalEqualities(x, x + size_.variables, c + size_.equalities);
      scatterEqualities(c);
      return 0;
    }

    /**
     * @brief A component owns its inequality rows, so the kernel writes the view
     */
    template <typename scalar_type, typename index_type>
    int Component<scalar_type, index_type>::evaluateInequalities()
    {
      ScalarT* x = local_x_.getData();
      gather(x);
      inequalities(x, x + size_.variables, d_.getData());
      return 0;
    }

    /// Local variables: internal, then those each port reads
    template <typename scalar_type, typename index_type>
    void Component<scalar_type, index_type>::gather(ScalarT* x) const
    {
      x = std::copy_n(x_.getData(), size_.variables, x);
      for (const Port& port : ports_)
      {
        x = std::copy_n(port.component->x_.getData() + port.offset.variables, port.size.variables, x);
      }
    }

    /// Global column of each local variable
    template <typename scalar_type, typename index_type>
    void Component<scalar_type, index_type>::variableIndices(IdxT* columns) const
    {
      std::iota(columns, columns + size_.variables, offset_.variables);
      columns += size_.variables;
      for (const Port& port : ports_)
      {
        std::iota(columns, columns + port.size.variables, port.component->offset_.variables + port.offset.variables);
        columns += port.size.variables;
      }
    }

    /// Global row of each local equality contribution
    template <typename scalar_type, typename index_type>
    void Component<scalar_type, index_type>::equalityIndices(IdxT* rows) const
    {
      std::iota(rows, rows + size_.equalities, offset_.equalities);
      rows += size_.equalities;
      for (const Port& port : ports_)
      {
        std::iota(rows, rows + port.size.equalities, port.component->offset_.equalities + port.offset.equalities);
        rows += port.size.equalities;
      }
    }

    /// Add the local gradient to the gradient of each variable it belongs to
    template <typename scalar_type, typename index_type>
    void Component<scalar_type, index_type>::scatterGradient(const RealT* gradient)
    {
      RealT* own = gradient_.getData();
      for (IdxT i = 0; i < size_.variables; ++i)
      {
        own[i] += *gradient++;
      }
      for (const Port& port : ports_)
      {
        RealT* target = port.component->gradient_.getData() + port.offset.variables;
        for (IdxT i = 0; i < port.size.variables; ++i)
        {
          target[i] += *gradient++;
        }
      }
    }

    /// Add the local equality contributions to the row each belongs to
    template <typename scalar_type, typename index_type>
    void Component<scalar_type, index_type>::scatterEqualities(const ScalarT* c)
    {
      ScalarT* own = c_.getData();
      for (IdxT i = 0; i < size_.equalities; ++i)
      {
        own[i] += *c++;
      }
      for (const Port& port : ports_)
      {
        ScalarT* target = port.component->c_.getData() + port.offset.equalities;
        for (IdxT i = 0; i < port.size.equalities; ++i)
        {
          target[i] += *c++;
        }
      }
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
