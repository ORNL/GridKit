/**
 * @file SystemModelDependencyTracking.cpp
 * @brief Optimal power flow system derivatives from `DependencyTracking::Variable`.
 *
 * @note Used for testing the Enzyme derivatives.
 */

#include "SystemModelImpl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace
    {
      /**
       * @brief Jacobian of `functions` from their dependencies
       *
       * @param[in] functions - Constraint functions evaluated at the current variables
       * @param[in] cols - Number of variables
       */
      template <typename CsrMatrixT, typename RealT, typename VectorT, typename IdxT>
      std::unique_ptr<CsrMatrixT> dependencyJacobian(const VectorT& functions, IdxT cols)
      {
        const IdxT  rows = functions.getSize();
        const auto* g    = functions.getData();

        IdxT* row_ptrs = new IdxT[rows + 1];
        row_ptrs[0]    = 0;
        for (IdxT row = 0; row < rows; ++row)
        {
          row_ptrs[row + 1] = row_ptrs[row] + g[row].getDependencies().size();
        }

        IdxT*  col_data = new IdxT[row_ptrs[rows]];
        RealT* vals     = new RealT[row_ptrs[rows]];

        IdxT counter = 0;
        for (IdxT row = 0; row < rows; ++row)
        {
          for (const auto& [col, value] : g[row].getDependencies())
          {
            col_data[counter] = col;
            vals[counter]     = value;
            ++counter;
          }
        }

        return std::make_unique<CsrMatrixT>(rows, cols, row_ptrs[rows], &row_ptrs, &col_data, &vals);
      }
    } // namespace

    /**
     * @brief Derivatives are built from dependencies on each evaluation
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::allocateDerivatives()
    {
      const int status = evaluateEqualityJacobian();
      if (status != 0)
      {
        return status;
      }
      return evaluateInequalityJacobian();
    }

    /**
     * @brief Gradient from the objective dependencies
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateGradient()
    {
      const int status = evaluateObjective();
      if (status != 0)
      {
        return status;
      }
      gradient_.setToZero();
      RealT* gradient = gradient_.getData();
      for (const auto& [variable, derivative] : f_.getDependencies())
      {
        gradient[variable] = derivative;
      }
      gradient_.setDataUpdated();

      return 0;
    }

    /**
     * @brief Equality constraint Jacobian from the dependencies
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateEqualityJacobian()
    {
      const int status = evaluateEqualities();
      if (status != 0)
      {
        return status;
      }
      equality_jacobian_ = dependencyJacobian<CsrMatrixT, RealT>(c_, x_.getSize());

      return 0;
    }

    /**
     * @brief Inequality constraint Jacobian from the dependencies
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateInequalityJacobian()
    {
      const int status = evaluateInequalities();
      if (status != 0)
      {
        return status;
      }
      inequality_jacobian_ = dependencyJacobian<CsrMatrixT, RealT>(d_, x_.getSize());

      return 0;
    }

    /**
     * @brief Not available, because dependencies hold first derivatives only
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateHessian([[maybe_unused]] RealT        sigma,
                                                              [[maybe_unused]] const RealT* lambda_c,
                                                              [[maybe_unused]] const RealT* lambda_d)
    {
      Log::error() << "OptimalPowerFlow::SystemModel: no Hessian for DependencyTracking::Variable\n";
      return 1;
    }

    // Available template instantiations
    template class SystemModel<DependencyTracking::Variable, size_t>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
