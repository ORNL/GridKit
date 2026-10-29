/**
 * @file ComponentImpl.hpp
 * @brief Definition of `Component`, included by model translation units.
 *
 * A component evaluates on a local copy. `gather` copies its own variables,
 * then the variables at each port, into one fixed-size array, the model
 * evaluates there, and the scatter functions add the results back into the
 * system vectors.
 *
 * Gathers only read, so they are safe in parallel. Scatters are not: every
 * device at a bus adds to that bus's entries, so those additions must be
 * synchronized, and buses with many devices are the likely bottleneck.
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
          || variable_indices_.resize(local.variables) != 0
          || constraint_indices_.resize(local.constraints) != 0
          || jacobian_entries_.resize(local.constraints * local.variables) != 0
          || hessian_entries_.resize(local.variables * local.variables) != 0)
      {
        Log::error() << id_ << ": failed to allocate the local buffers\n";
        throw std::bad_alloc{};
      }
    }

    template <typename scalar_type, typename index_type>
    int Component<scalar_type, index_type>::connect(IdxT port, Component& component, const Block& offset)
    {
      if (bound_ || port >= ports_.size())
      {
        return 1;
      }

      const Block& size      = ports_[port].size;
      const Block& available = component.size_;
      if (offset.variables > available.variables
          || size.variables > available.variables - offset.variables
          || offset.constraints > available.constraints
          || size.constraints > available.constraints - offset.constraints)
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

    /// Copy the internal variables, then the variables each port reads from its component, into `x`
    template <typename scalar_type, typename index_type>
    void Component<scalar_type, index_type>::gather(ScalarT* x) const
    {
      x = std::copy_n(x_.getData(), size_.variables, x);
      for (const Port& port : ports_)
      {
        x = std::copy_n(port.component->x_.getData() + port.offset.variables, port.size.variables, x);
      }
    }

    /// Multipliers of the internal, then each port's constraints
    template <typename scalar_type, typename index_type>
    void Component<scalar_type, index_type>::gatherMultipliers(RealT* lambda) const
    {
      lambda = std::copy_n(lambda_.getData(), size_.constraints, lambda);
      for (const Port& port : ports_)
      {
        lambda = std::copy_n(port.component->lambda_.getData() + port.offset.constraints, port.size.constraints, lambda);
      }
    }

    /// Global column of each local variable and row of each local constraint contribution
    template <typename scalar_type, typename index_type>
    void Component<scalar_type, index_type>::setIndices()
    {
      IdxT* columns = variable_indices_.getData();
      IdxT* rows    = constraint_indices_.getData();
      std::iota(columns, columns + size_.variables, offset_.variables);
      std::iota(rows, rows + size_.constraints, offset_.constraints);
      columns += size_.variables;
      rows    += size_.constraints;
      for (const Port& port : ports_)
      {
        std::iota(columns, columns + port.size.variables, port.component->offset_.variables + port.offset.variables);
        std::iota(rows, rows + port.size.constraints, port.component->offset_.constraints + port.offset.constraints);
        columns += port.size.variables;
        rows    += port.size.constraints;
      }
    }

    /// Add the local gradient to the system gradient
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

    /// Add the local contributions to the system constraints
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
