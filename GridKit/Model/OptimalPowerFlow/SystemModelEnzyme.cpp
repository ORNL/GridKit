/**
 * @file SystemModelEnzyme.cpp
 * @brief Optimal power flow system derivatives from component Enzyme entries.
 */

#include <algorithm>
#include <new>

#include <GridKit/LinearAlgebra/SparseMatrix/CooMatrix.hpp>

#include "SystemModelImpl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateGradient()
    {
      gradient_.setToZero();
      for (auto& component : components_)
      {
        const int status = component->evaluateGradient();
        if (status != 0)
        {
          return status;
        }
      }
      gradient_.setDataUpdated();

      return 0;
    }

    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateJacobian()
    {
      for (auto& component : components_)
      {
        const int status = component->evaluateJacobian();
        if (status != 0)
        {
          return status;
        }
      }

      return refill(*jacobian_, jacobian_assembly_);
    }

    /**
     * @brief Lower triangle of the Hessian of \f$\sigma f + \lambda^T g\f$ at
     * the multipliers `lambda()`
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateHessian(RealT sigma)
    {
      for (auto& component : components_)
      {
        const int status = component->evaluateHessian(sigma);
        if (status != 0)
        {
          return status;
        }
      }

      return refill(*hessian_, hessian_assembly_);
    }

    /**
     * @brief Construct the CSR patterns from the component entries
     *
     * Component entries are structural, so the patterns hold at every point
     * and every multiplier, including zero.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::allocateDerivatives()
    {
      jacobian_assembly_.entries.clear();
      hessian_assembly_.entries.clear();
      for (auto& component : components_)
      {
        if (component->evaluateJacobian() != 0 || component->evaluateHessian(ONE<RealT>) != 0)
        {
          return 1;
        }

        jacobian_assembly_.entries.push_back(&component->jacobianEntries());
        hessian_assembly_.entries.push_back(&component->hessianEntries());
      }

      jacobian_ = assemble(g_.getSize(), x_.getSize(), jacobian_assembly_);
      hessian_  = assemble(x_.getSize(), x_.getSize(), hessian_assembly_);

      return 0;
    }

    /**
     * @brief Sort and deduplicate the component entries into a CSR matrix
     *
     * @param[in] rows - Number of rows
     * @param[in] cols - Number of columns
     * @param[in,out] assembly - Component entries, and the CSR entry of each on return
     */
    template <typename scalar_type, typename index_type>
    std::unique_ptr<typename SystemModel<scalar_type, index_type>::CsrMatrixT>
    SystemModel<scalar_type, index_type>::assemble(IdxT rows, IdxT cols, Assembly& assembly) const
    {
      IdxT nnz_dup = 0;
      for (const CooBufferT* entries : assembly.entries)
      {
        nnz_dup += entries->nnz;
      }

      if (assembly.map_to_csr.resize(nnz_dup) != 0)
      {
        throw std::bad_alloc{};
      }

      // COO triplet arrays, handed off to the COO matrix
      IdxT*  rows_dup = new IdxT[nnz_dup];
      IdxT*  cols_dup = new IdxT[nnz_dup];
      RealT* vals_dup = new RealT[nnz_dup];

      IdxT counter = 0;
      for (const CooBufferT* entries : assembly.entries)
      {
        std::copy_n(entries->rows.getData(), entries->nnz, rows_dup + counter);
        std::copy_n(entries->cols.getData(), entries->nnz, cols_dup + counter);
        std::copy_n(entries->values.getData(), entries->nnz, vals_dup + counter);
        counter += entries->nnz;
      }

      LinearAlgebra::CooMatrix<RealT, IdxT> coo(rows, cols, nnz_dup, &rows_dup, &cols_dup, &vals_dup);

      // Populate CSR data with sort and deduplicate
      IdxT*      row_ptrs = coo.getCsrRowData();
      const IdxT nnz      = coo.getNnz();

      IdxT*  csr_cols = new IdxT[nnz];
      RealT* csr_vals = new RealT[nnz];
      std::copy(coo.getColData(), coo.getColData() + nnz, csr_cols);
      std::copy(coo.getValues(), coo.getValues() + nnz, csr_vals);

      const IdxT* map_to_sorted = coo.getMapToSorted();
      const IdxT* map_to_dedup  = coo.getMapToDeduplicated();

      IdxT* map = assembly.map_to_csr.getData();
      for (IdxT i = 0; i < nnz_dup; ++i)
      {
        map[map_to_sorted[i]] = map_to_dedup[i];
      }
      assembly.map_to_csr.setDataUpdated();

      return std::make_unique<CsrMatrixT>(rows, cols, nnz, &row_ptrs, &csr_cols, &csr_vals);
    }

    /**
     * @brief Sum the component entries into the CSR values
     *
     * @return 1 if the number of component entries changed since assembly
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::refill(CsrMatrixT& matrix, const Assembly& assembly) const
    {
      IdxT nnz_dup = 0;
      for (const CooBufferT* entries : assembly.entries)
      {
        nnz_dup += entries->nnz;
      }

      if (nnz_dup != assembly.map_to_csr.getSize())
      {
        Log::error() << "OptimalPowerFlow::SystemModel: " << nnz_dup << " derivative entries, but "
                     << assembly.map_to_csr.getSize() << " at assembly\n";
        return 1;
      }

      const IdxT* map    = assembly.map_to_csr.getData();
      RealT*      values = matrix.getValues();
      std::fill_n(values, matrix.getNnz(), ZERO<RealT>);

      IdxT counter = 0;
      for (const CooBufferT* entries : assembly.entries)
      {
        const RealT* entry_values = entries->values.getData();
        for (IdxT k = 0; k < entries->nnz; ++k)
        {
          values[map[counter++]] += entry_values[k];
        }
      }

      return matrix.setUpdated(memory::HOST);
    }

    // Available template instantiations
    template class SystemModel<double, size_t>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
