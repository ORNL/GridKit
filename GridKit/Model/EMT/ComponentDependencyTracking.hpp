/**
 * @file ComponentDependencyTracking.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Partials of EMT models whose parameters size their layouts.
 *
 * The value and derivative of local value j are tracked as 2j and 2j + 1, at
 * their current values. Dependencies scaled by zero are kept, so the pattern
 * does not change. An output that is not connected is not tracked and has no
 * entries.
 *
 * A linear model forwards `tagDifferentiable`, which stores its partials. A
 * model that is not linear also forwards `evaluateJacobian`, which stores
 * them again at every Jacobian evaluation, without allocating the entries
 * again.
 */

#pragma once

#include <algorithm>
#include <new>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Store the entries and the partials dF/dy and dF/dy' of `model` at
     *        the current values.
     *
     * @tparam ModelT Model providing the kernels.
     * @param[in] model This component as its derived type.
     * @return 0.
     *
     * @throws std::bad_alloc if the entries cannot be allocated.
     *
     * @pre The system has bound this component and resolved its indices.
     */
    template <typename scalar_type, typename index_type>
    template <typename ModelT>
    int Component<scalar_type, index_type>::tagDifferentiable(const ModelT& model)
    {
      using Variable = DependencyTracking::Variable;

      const IdxT n = variable_indices_.getSize();
      const IdxT m = residual_indices_.size();

      // A model that is not linear tracks at every Jacobian evaluation, so the components of a model share these
      static thread_local std::vector<Variable> y;
      static thread_local std::vector<Variable> yp;
      static thread_local std::vector<Variable> f;

      const bool numbered = y.size() == n;
      if (!numbered)
      {
        y.assign(n, Variable{});
        yp.assign(n, Variable{});
        for (IdxT j = 0; j < n; ++j)
        {
          y[j].setVariableNumber(2 * j);
          yp[j].setVariableNumber(2 * j + 1);
        }
      }
      f.assign(m, Variable{});

      const auto [w, wp] = gather();
      for (IdxT j = 0; j < size_; ++j)
      {
        y[j].setValue(static_cast<double>(y_.getData()[j]));
        yp[j].setValue(static_cast<double>(yp_.getData()[j]));
      }
      for (IdxT j = size_; j < n; ++j)
      {
        y[j].setValue(static_cast<double>(w[j - size_]));
        yp[j].setValue(static_cast<double>(wp[j - size_]));
      }
      evaluateKernels(model, y.data(), yp.data(), y.data() + size_, yp.data() + size_, f.data());

      // The value and derivative of a local variable are adjacent in the dependency map
      const auto partials = [this, m](auto entry)
      {
        for (IdxT row = 0; row < m; ++row)
        {
          if (residual_indices_[row] == INVALID_INDEX<IdxT>)
          {
            continue;
          }

          IdxT  last  = INVALID_INDEX<IdxT>;
          RealT dfdy  = ZERO<RealT>;
          RealT dfdyp = ZERO<RealT>;
          for (const auto& [number, value] : f[row].getDependencies())
          {
            const auto local = static_cast<IdxT>(number / 2);
            if (local != last && last != INVALID_INDEX<IdxT>)
            {
              entry(row, last, dfdy, dfdyp);
              dfdy  = ZERO<RealT>;
              dfdyp = ZERO<RealT>;
            }
            last = local;

            if (number % 2 == 0)
            {
              dfdy = static_cast<RealT>(value);
            }
            else
            {
              dfdyp = static_cast<RealT>(value);
            }
          }
          if (last != INVALID_INDEX<IdxT>)
          {
            entry(row, last, dfdy, dfdyp);
          }
        }
      };

      // The pattern does not change, so the entries are allocated once
      IdxT nnz = 0;
      partials([&nnz](IdxT, IdxT, RealT, RealT)
               { ++nnz; });
      if (jacobian_entries_.rows.getSize() != nnz || dfdy_.getSize() != nnz || dfdyp_.getSize() != nnz)
      {
        if (jacobian_entries_.resize(nnz) != 0 || dfdy_.resize(nnz) != 0 || dfdyp_.resize(nnz) != 0)
        {
          Log::error() << id_ << ": failed to allocate the derivative entries\n";
          throw std::bad_alloc{};
        }
      }

      const IdxT* indices = variable_indices_.getData();
      IdxT*       rows    = jacobian_entries_.rows.getData();
      IdxT*       cols    = jacobian_entries_.cols.getData();
      RealT*      dfdy    = dfdy_.getData();
      RealT*      dfdyp   = dfdyp_.getData();
      IdxT        k       = 0;
      partials([&](IdxT row, IdxT local, RealT value, RealT derivative)
               {
                 rows[k]  = residual_indices_[row];
                 cols[k]  = indices[local];
                 dfdy[k]  = value;
                 dfdyp[k] = derivative;
                 ++k; });
      jacobian_entries_.nnz = nnz;
      return 0;
    }

    /**
     * @brief Store the partials of `model` at the current values, then set
     *        the Jacobian entries to dF/dy + alpha dF/dy'.
     *
     * A model that is not linear forwards `evaluateJacobian` here.
     *
     * @tparam ModelT Model providing the kernels.
     * @param[in] model This component as its derived type.
     * @return 0 on success, or the status of `tagDifferentiable`.
     *
     * @pre The system has bound this component and resolved its indices.
     */
    template <typename scalar_type, typename index_type>
    template <typename ModelT>
    int Component<scalar_type, index_type>::evaluateJacobian(const ModelT& model)
    {
      const int status = tagDifferentiable(model);
      if (status != 0)
      {
        return status;
      }
      // The refill of this class, not the override that forwards here
      return Component::evaluateJacobian();
    }
  } // namespace EMT
} // namespace GridKit
