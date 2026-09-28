/**
 * @file SystemModelEnzyme.cpp
 * @brief Optimal power flow system derivatives from component Enzyme entries.
 */

#include <algorithm>

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
      const ScalarT*           x = x_.getData();
      const std::vector<RealT> lambda(g_.getSize(), ZERO<RealT>);
      for (auto& component : components_)
      {
        component->evaluateJacobian(x);
        component->evaluateHessian(x, ONE<RealT>, lambda.data());
      }

      jacobian_ = assemble(g_.getSize(), x_.getSize(), &ComponentT::jacobian, jacobian_map_);
      hessian_  = assemble(x_.getSize(), x_.getSize(), &ComponentT::hessian, hessian_map_);

      return 0;
    }

    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateGradient()
    {
      const ScalarT* x = x_.getData();

      gradient_.setToZero();
      RealT* gradient = gradient_.getData();
      for (auto& component : components_)
      {
        component->evaluateGradient(x, gradient);
      }
      gradient_.setDataUpdated();

      return 0;
    }

    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateJacobian()
    {
      const ScalarT* x = x_.getData();
      for (auto& component : components_)
      {
        component->evaluateJacobian(x);
      }

      return refill(*jacobian_, &ComponentT::jacobian, jacobian_map_);
    }

    /**
     * @brief Lower triangle of the Hessian of \f$\sigma f + \lambda^T g\f$
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateHessian(RealT sigma, const RealT* lambda)
    {
      const ScalarT* x = x_.getData();
      for (auto& component : components_)
      {
        component->evaluateHessian(x, sigma, lambda);
      }

      return refill(*hessian_, &ComponentT::hessian, hessian_map_);
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
    SystemModel<scalar_type, index_type>::assemble(IdxT               rows,
                                                   IdxT               cols,
                                                   EntriesT           entries,
                                                   std::vector<IdxT>& map_to_csr) const
    {
      IdxT nnz_dup = 0;
      for (const auto& component : components_)
      {
        nnz_dup += ((*component).*entries)().values.size();
      }

      // COO triplet arrays, handed off to the COO matrix
      IdxT*  rows_dup = new IdxT[nnz_dup];
      IdxT*  cols_dup = new IdxT[nnz_dup];
      RealT* vals_dup = new RealT[nnz_dup];

      IdxT counter = 0;
      for (const auto& component : components_)
      {
        const CooEntriesT& coo = ((*component).*entries)();
        std::copy(coo.rows.begin(), coo.rows.end(), rows_dup + counter);
        std::copy(coo.cols.begin(), coo.cols.end(), cols_dup + counter);
        std::copy(coo.values.begin(), coo.values.end(), vals_dup + counter);
        counter += coo.values.size();
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

      map_to_csr.resize(nnz_dup);
      for (IdxT i = 0; i < nnz_dup; ++i)
      {
        map_to_csr[map_to_sorted[i]] = map_to_dedup[i];
      }

      return std::make_unique<CsrMatrixT>(rows, cols, nnz, &row_ptrs, &csr_cols, &csr_vals);
    }

    /**
     * @brief Sum the component entries into the CSR values
     *
     * @return 1 if the number of component entries changed since assembly
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::refill(CsrMatrixT&              matrix,
                                                     EntriesT                 entries,
                                                     const std::vector<IdxT>& map_to_csr) const
    {
      IdxT nnz_dup = 0;
      for (const auto& component : components_)
      {
        nnz_dup += ((*component).*entries)().values.size();
      }

      if (nnz_dup != map_to_csr.size())
      {
        Log::error() << "OptimalPowerFlow::SystemModel: " << nnz_dup << " derivative entries, but "
                     << map_to_csr.size() << " at assembly\n";
        return 1;
      }

      RealT* values = matrix.getValues();
      std::fill_n(values, matrix.getNnz(), ZERO<RealT>);

      IdxT counter = 0;
      for (const auto& component : components_)
      {
        for (const RealT value : ((*component).*entries)().values)
        {
          values[map_to_csr[counter]] += value;
          ++counter;
        }
      }

      return 0;
    }

    // Available template instantiations
    template class SystemModel<double, size_t>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
