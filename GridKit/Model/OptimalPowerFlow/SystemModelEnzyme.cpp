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
    /**
     * @brief Construct the CSR patterns from the component entries at the
     * starting point
     *
     * Component entries are structural, so the patterns hold at every point
     * and every multiplier, including zero.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::allocateDerivatives()
    {
      RealVectorT lambda_c;
      RealVectorT lambda_d;
      if (lambda_c.resize(c_.getSize()) != 0 || lambda_d.resize(d_.getSize()) != 0)
      {
        return 1;
      }
      for (auto& component : components_)
      {
        int status = component->evaluateEqualityJacobian();
        if (status != 0)
        {
          return status;
        }

        status = component->evaluateInequalityJacobian();
        if (status != 0)
        {
          return status;
        }

        status = component->evaluateHessian(ONE<RealT>, lambda_c.getData(), lambda_d.getData());
        if (status != 0)
        {
          return status;
        }
      }

      equality_jacobian_   = assemble(c_.getSize(), x_.getSize(), &ComponentT::equalityEntries, equality_jacobian_map_);
      inequality_jacobian_ = assemble(d_.getSize(), x_.getSize(), &ComponentT::inequalityEntries, inequality_jacobian_map_);
      hessian_             = assemble(x_.getSize(), x_.getSize(), &ComponentT::hessianEntries, hessian_map_);

      return 0;
    }

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
    int SystemModel<scalar_type, index_type>::evaluateEqualityJacobian()
    {
      for (auto& component : components_)
      {
        const int status = component->evaluateEqualityJacobian();
        if (status != 0)
        {
          return status;
        }
      }

      return refill(*equality_jacobian_, &ComponentT::equalityEntries, equality_jacobian_map_);
    }

    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateInequalityJacobian()
    {
      for (auto& component : components_)
      {
        const int status = component->evaluateInequalityJacobian();
        if (status != 0)
        {
          return status;
        }
      }

      return refill(*inequality_jacobian_, &ComponentT::inequalityEntries, inequality_jacobian_map_);
    }

    /**
     * @brief Lower triangle of the Hessian of
     * \f$\sigma f + \lambda_c^T c + \lambda_d^T d\f$
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateHessian(RealT sigma, const RealT* lambda_c, const RealT* lambda_d)
    {
      for (auto& component : components_)
      {
        const int status = component->evaluateHessian(sigma, lambda_c, lambda_d);
        if (status != 0)
        {
          return status;
        }
      }

      return refill(*hessian_, &ComponentT::hessianEntries, hessian_map_);
    }

    /**
     * @brief Sort and deduplicate the component entries into a CSR matrix
     *
     * @param[in] rows - Number of rows
     * @param[in] cols - Number of columns
     * @param[in] entries - Component entries to assemble
     * @param[out] map_to_csr - CSR index of each component entry
     */
    template <typename scalar_type, typename index_type>
    std::unique_ptr<typename SystemModel<scalar_type, index_type>::CsrMatrixT>
    SystemModel<scalar_type, index_type>::assemble(IdxT          rows,
                                                   IdxT          cols,
                                                   EntriesT      entries,
                                                   IndexVectorT& map_to_csr) const
    {
      IdxT nnz_dup = 0;
      for (const auto& component : components_)
      {
        nnz_dup += ((*component).*entries)().nnz;
      }

      if (map_to_csr.resize(nnz_dup) != 0)
      {
        throw std::bad_alloc{};
      }

      // COO triplet arrays, handed off to the COO matrix
      IdxT*  rows_dup = new IdxT[nnz_dup];
      IdxT*  cols_dup = new IdxT[nnz_dup];
      RealT* vals_dup = new RealT[nnz_dup];

      IdxT counter = 0;
      for (const auto& component : components_)
      {
        const CooBufferT& buffer = ((*component).*entries)();
        std::copy_n(buffer.rows.getData(), buffer.nnz, rows_dup + counter);
        std::copy_n(buffer.cols.getData(), buffer.nnz, cols_dup + counter);
        std::copy_n(buffer.values.getData(), buffer.nnz, vals_dup + counter);
        counter += buffer.nnz;
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

      IdxT* map = map_to_csr.getData();
      for (IdxT i = 0; i < nnz_dup; ++i)
      {
        map[map_to_sorted[i]] = map_to_dedup[i];
      }
      map_to_csr.setDataUpdated();

      return std::make_unique<CsrMatrixT>(rows, cols, nnz, &row_ptrs, &csr_cols, &csr_vals);
    }

    /**
     * @brief Sum the component entries into the CSR values
     *
     * @return 1 if the number of component entries changed since assembly
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::refill(CsrMatrixT&         matrix,
                                                     EntriesT            entries,
                                                     const IndexVectorT& map_to_csr) const
    {
      IdxT nnz_dup = 0;
      for (const auto& component : components_)
      {
        nnz_dup += ((*component).*entries)().nnz;
      }

      if (nnz_dup != map_to_csr.getSize())
      {
        Log::error() << "OptimalPowerFlow::SystemModel: " << nnz_dup << " derivative entries, but "
                     << map_to_csr.getSize() << " at assembly\n";
        return 1;
      }

      const IdxT* map    = map_to_csr.getData();
      RealT*      values = matrix.getValues();
      std::fill_n(values, matrix.getNnz(), ZERO<RealT>);

      IdxT counter = 0;
      for (const auto& component : components_)
      {
        const CooBufferT& buffer           = ((*component).*entries)();
        const RealT*      component_values = buffer.values.getData();
        for (IdxT k = 0; k < buffer.nnz; ++k)
        {
          values[map[counter++]] += component_values[k];
        }
      }

      return matrix.setUpdated(memory::HOST);
    }

    // Available template instantiations
    template class SystemModel<double, size_t>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
