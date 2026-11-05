/**
 * @file SystemModel.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief EMT system derivatives assembled from the component entries.
 */

#include <algorithm>
#include <new>
#include <vector>

#include <GridKit/LinearAlgebra/SparseMatrix/CooMatrix.hpp>

#include "SystemModelImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateJacobian()
    {
      for (ComponentT* component : blocks_)
      {
        const int status = component->evaluateJacobian();
        if (status != 0)
        {
          return status;
        }
      }

      return refill(*jacobian_, jacobian_assembly_);
    }

    /// A variable is differential when a residual depends on its derivative
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::tagDifferentiable()
    {
      tag_.assign(y_.getSize(), false);
      for (ComponentT* component : blocks_)
      {
        const int status = component->tagDifferentiable();
        if (status != 0)
        {
          return status;
        }

        const IdxT*  cols  = component->jacobian_entries_.cols.getData();
        const RealT* dfdyp = component->dfdyp_.getData();
        for (IdxT k = 0; k < component->dfdyp_.getSize(); ++k)
        {
          if (dfdyp[k] != ZERO<RealT>)
          {
            tag_[cols[k]] = true;
          }
        }
      }
      return 0;
    }

    /// CSR pattern from the component entries, which are structural
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::allocateDerivatives()
    {
      jacobian_assembly_.entries.clear();
      for (ComponentT* component : blocks_)
      {
        if (component->tagDifferentiable() != 0 || component->evaluateJacobian() != 0)
        {
          return 1;
        }
        jacobian_assembly_.entries.push_back(&component->jacobian_entries_);
      }

      jacobian_ = assemble(y_.getSize(), y_.getSize(), jacobian_assembly_);
      return 0;
    }

    /// Sort and deduplicate the entries into a CSR matrix, and map each entry to its CSR value
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
        Log::error() << "EMT::SystemModel: failed to allocate the CSR map\n";
        throw std::bad_alloc{};
      }

      // Handed off to the COO matrix
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

    /// Sum the component entries into the CSR values
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
        Log::error() << "EMT::SystemModel: " << nnz_dup << " derivative entries, but "
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
    template class Component<double, size_t>;
    template class SystemModel<double, size_t>;
  } // namespace EMT
} // namespace GridKit
