/**
 * @file ComponentDependencyTracking.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Partials of EMT models that are linear between events.
 *
 * The value and derivative of local variable j are tracked as 2j and 2j + 1.
 * Dependencies scaled by zero are kept, so events keep the pattern.
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
    template <typename scalar_type, typename index_type>
    template <typename ModelT>
    int Component<scalar_type, index_type>::tagDifferentiable(const ModelT& model)
    {
      using Variable = DependencyTracking::Variable;

      const IdxT            n = variable_indices_.getSize();
      std::vector<Variable> y(n);
      std::vector<Variable> yp(n);
      std::vector<Variable> f(n);
      for (IdxT j = 0; j < n; ++j)
      {
        y[j].setVariableNumber(2 * j);
        yp[j].setVariableNumber(2 * j + 1);
      }
      model.internalResidual(y.data(), yp.data(), y.data() + size_, yp.data() + size_, f.data());
      model.externalResidual(y.data(), yp.data(), y.data() + size_, yp.data() + size_, f.data() + size_);

      const IdxT*        indices = variable_indices_.getData();
      std::vector<IdxT>  rows;
      std::vector<IdxT>  cols;
      std::vector<RealT> dfdy;
      std::vector<RealT> dfdyp;
      for (IdxT row = 0; row < n; ++row)
      {
        // The value and derivative of a local variable are adjacent in the dependency map
        IdxT last = INVALID_INDEX<IdxT>;
        for (const auto& [number, value] : f[row].getDependencies())
        {
          const auto local = static_cast<IdxT>(number / 2);
          if (local != last)
          {
            rows.push_back(indices[row]);
            cols.push_back(indices[local]);
            dfdy.push_back(ZERO<RealT>);
            dfdyp.push_back(ZERO<RealT>);
            last = local;
          }

          if (number % 2 == 0)
          {
            dfdy.back() = static_cast<RealT>(value);
          }
          else
          {
            dfdyp.back() = static_cast<RealT>(value);
          }
        }
      }

      const auto nnz = static_cast<IdxT>(rows.size());
      if (jacobian_entries_.resize(nnz) != 0 || dfdy_.resize(nnz) != 0 || dfdyp_.resize(nnz) != 0)
      {
        Log::error() << id_ << ": failed to allocate the derivative entries\n";
        throw std::bad_alloc{};
      }
      std::copy(rows.begin(), rows.end(), jacobian_entries_.rows.getData());
      std::copy(cols.begin(), cols.end(), jacobian_entries_.cols.getData());
      std::copy(dfdy.begin(), dfdy.end(), dfdy_.getData());
      std::copy(dfdyp.begin(), dfdyp.end(), dfdyp_.getData());
      jacobian_entries_.nnz = nnz;
      return 0;
    }
  } // namespace EMT
} // namespace GridKit
