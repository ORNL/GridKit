#include <algorithm>
#include <cmath>
#include <cstddef>

#include <GridKit/AutomaticDifferentiation/Enzyme/Optimization/ConstraintJacobian.hpp>
#include <GridKit/AutomaticDifferentiation/Enzyme/Optimization/LagrangianHessian.hpp>

#include "BranchImpl.hpp"

namespace GridKit
{
  namespace OPF
  {
    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::evaluateJacobian()
    {
      const int status = GridKit::Enzyme::Optimization::
          ConstraintJacobian<Branch, IdxT>::eval(
              this,
              VARIABLE_COUNT,
              jacobianEntries(),
              variables_.data(),
              jacobian_values_.data());
      if (status != 0)
      {
        return status;
      }
      for (const RealT value : jacobianValues())
      {
        if (!std::isfinite(value))
        {
          return 1;
        }
      }
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::evaluateHessian(
        RealT                  objective_factor,
        std::span<const RealT> local_multipliers)
    {
      if (local_multipliers.size() != constraint_indices_.size())
      {
        return 1;
      }
      local_multipliers_.fill(RealT{0});
      std::copy(local_multipliers.begin(),
                local_multipliers.end(),
                local_multipliers_.begin());

      const int status = GridKit::Enzyme::Optimization::
          LagrangianHessian<Branch, IdxT>::eval(
              this,
              VARIABLE_COUNT,
              hessianEntries(),
              variables_.data(),
              objective_factor,
              local_multipliers_.data(),
              hessian_values_.data());
      if (status != 0)
      {
        return status;
      }
      for (const RealT value : hessian_values_)
      {
        if (!std::isfinite(value))
        {
          return 1;
        }
      }
      return 0;
    }

    template class Branch<double, int>;
    template class Branch<double, long int>;
    template class Branch<double, std::size_t>;
  } // namespace OPF
} // namespace GridKit
