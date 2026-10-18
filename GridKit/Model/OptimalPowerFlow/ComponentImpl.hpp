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
      tag_.assign(size_.constraints, false);
      for (IdxT port = 0; port < ModelT::NUM_PORT_VARIABLES.size(); ++port)
      {
        ports_.push_back({{ModelT::NUM_PORT_VARIABLES[port], ModelT::NUM_PORT_CONSTRAINTS[port]}});
      }

      constexpr Block local = localSize<ModelT>();
      if (local_x_.resize(local.variables) != 0
          || local_g_.resize(local.constraints) != 0
          || jacobian_entries_.resize(local.constraints * local.variables) != 0
          || hessian_entries_.resize(local.variables * local.variables) != 0)
      {
        throw std::bad_alloc{};
      }
    }

    /**
     * @brief Read the variables and contribute to the constraints of
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
          || !fits(offset.constraints, size.constraints, component.size_.constraints))
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
    int Component<scalar_type, index_type>::evaluateConstraints()
    {
      ScalarT* x = local_x_.getData();
      ScalarT* g = local_g_.getData();
      gather(x);
      local_g_.setToZero();
      internalConstraints(x, x + size_.variables, g);
      externalConstraints(x, x + size_.variables, g + size_.constraints);
      scatterConstraints(g);
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

    /// Multipliers of the local constraint contributions: internal, then those at each port
    template <typename scalar_type, typename index_type>
    void Component<scalar_type, index_type>::gatherMultipliers(RealT* lambda) const
    {
      lambda = std::copy_n(lambda_.getData(), size_.constraints, lambda);
      for (const Port& port : ports_)
      {
        lambda = std::copy_n(port.component->lambda_.getData() + port.offset.constraints, port.size.constraints, lambda);
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

    /// Global row of each local constraint contribution
    template <typename scalar_type, typename index_type>
    void Component<scalar_type, index_type>::constraintIndices(IdxT* rows) const
    {
      std::iota(rows, rows + size_.constraints, offset_.constraints);
      rows += size_.constraints;
      for (const Port& port : ports_)
      {
        std::iota(rows, rows + port.size.constraints, port.component->offset_.constraints + port.offset.constraints);
        rows += port.size.constraints;
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

    /// Add the local constraint contributions to the constraint each belongs to
    template <typename scalar_type, typename index_type>
    void Component<scalar_type, index_type>::scatterConstraints(const ScalarT* g)
    {
      ScalarT* own = g_.getData();
      for (IdxT i = 0; i < size_.constraints; ++i)
      {
        own[i] += *g++;
      }
      for (const Port& port : ports_)
      {
        ScalarT* target = port.component->g_.getData() + port.offset.constraints;
        for (IdxT i = 0; i < port.size.constraints; ++i)
        {
          target[i] += *g++;
        }
      }
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
