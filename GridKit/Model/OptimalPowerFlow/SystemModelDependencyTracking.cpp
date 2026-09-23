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
    /**
     * @brief Derivatives are built from dependencies on each evaluation
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::allocateDerivatives()
    {
      return 0;
    }

    /**
     * @brief Gradient from the objective dependencies
     *
     * @pre `evaluateObjective()` at the current variables
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateGradient()
    {
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
     * @brief Jacobian from the constraint dependencies
     *
     * @pre `evaluateConstraints()` at the current variables
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateJacobian()
    {
      const IdxT     m = g_.getSize();
      const ScalarT* g = g_.getData();

      IdxT* row_ptrs = new IdxT[m + 1];
      row_ptrs[0]    = 0;
      for (IdxT row = 0; row < m; ++row)
      {
        row_ptrs[row + 1] = row_ptrs[row] + g[row].getDependencies().size();
      }

      IdxT*  cols = new IdxT[row_ptrs[m]];
      RealT* vals = new RealT[row_ptrs[m]];

      IdxT counter = 0;
      for (IdxT row = 0; row < m; ++row)
      {
        for (const auto& [col, value] : g[row].getDependencies())
        {
          cols[counter] = col;
          vals[counter] = value;
          ++counter;
        }
      }

      jacobian_ = std::make_unique<CsrMatrixT>(m, x_.getSize(), row_ptrs[m], &row_ptrs, &cols, &vals);

      return 0;
    }

    /**
     * @brief Not available, because dependencies hold first derivatives only
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateHessian([[maybe_unused]] RealT        sigma,
                                                              [[maybe_unused]] const RealT* lambda)
    {
      Log::error() << "OptimalPowerFlow::SystemModel: no Hessian for DependencyTracking::Variable\n";
      return 1;
    }

    // Available template instantiations
    template class SystemModel<DependencyTracking::Variable, size_t>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
