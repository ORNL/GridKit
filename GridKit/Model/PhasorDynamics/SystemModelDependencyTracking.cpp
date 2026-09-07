#include "SystemModelImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief By default, Jacobians are not available
     *
     * DependencyTracking::Variable stores the Jacobian as dependency maps,
     * updated during calls to evaluateResidual().
     */
    template <typename scalar_type, typename index_type>
    bool SystemModel<scalar_type, index_type>::hasJacobian()
    {
      Log::warning() << "GridKit was not built with Enzyme. "
                     << "DependencyTracking::Variable Jacobians are only available for testing in PhasorDynamics.\n";

      return false;
    }

    /**
     * @brief Evaluate system DependencyTracking::Variable Jacobian.
     *
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateJacobian()
    {
      using DependencyMap = typename ScalarT::DependencyMap;
    
      const auto* f = f_.getData();
    
      if (csr_jac_ == nullptr)
      {
        IdxT* row_ptrs = new IdxT[static_cast<size_t>(size_) + 1];
        row_ptrs[0] = 0;
    
        // Count the number of non-zeros
        IdxT nnz = 0;
        for (IdxT row = 0; row < size_; ++row)
        {
          DependencyMap row_map;
    
          for (const auto& dep : f[row].getDependencies())
          {
            const auto col = dep.first;
    
            // Merge-count y and yp dependencies
            const IdxT jac_col = static_cast<IdxT>(col / 2);
    
            if (row_map.insert({jac_col, RealT{}}).second)
            {
              ++nnz;
            }
          }
    
          row_ptrs[static_cast<size_t>(row) + 1] = nnz;
        }
    
        // Allocate column and value pointers
        IdxT* cols  = new IdxT[static_cast<size_t>(nnz)];
        RealT* vals = new RealT[static_cast<size_t>(nnz)];
    
        // Store column and values
        IdxT i = 0;
        for (IdxT row = 0; row < size_; ++row)
        {
          DependencyMap row_map;
    
          for (const auto& dep : f[row].getDependencies())
          {
            const auto col = dep.first;
    
            const IdxT jac_col = static_cast<IdxT>(col / 2);
            // Even indices for y and odd indices for yp
            if (col % 2 == 0)
            {
              row_map[jac_col] += static_cast<RealT>(dep.second);
            }
            else
            {
              row_map[jac_col] += alpha_ * static_cast<RealT>(dep.second);
            }
          }
    
          for (const auto& entry : row_map)
          {
            cols[i] = static_cast<IdxT>(entry.first);
            vals[i] = static_cast<RealT>(entry.second);
            ++i;
          }
        }
    
        nnz_ = nnz;
        csr_jac_ = new CsrMatrixT(size_, size_, nnz_, &row_ptrs, &cols, &vals);
      }
      else
      {
        RealT* vals = csr_jac_->getValues();
    
        IdxT i = 0;
        for (IdxT row = 0; row < size_; ++row)
        {
          DependencyMap row_map;
    
          for (const auto& dep : f[row].getDependencies())
          {
            const auto col = dep.first;
    
            const IdxT jac_col = static_cast<IdxT>(col / 2);
            // Even indices for y and odd indices for yp
            if (col % 2 == 0)
            {
              row_map[jac_col] += static_cast<RealT>(dep.second);
            }
            else
            {
              row_map[jac_col] += alpha_ * static_cast<RealT>(dep.second);
            }
          }
    
          for (const auto& entry : row_map)
          {
            vals[i] = static_cast<RealT>(entry.second);
            ++i;
          }
        }
      }

      //Log::misc() << "System DependencyTracking Jacobian\n";
      //csr_jac_->print(Log::misc());
    
      return 0;
    }

    // Available template instantiations
    // template class SystemModel<DependencyTracking::Variable, long int>;
    template class SystemModel<DependencyTracking::Variable, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
