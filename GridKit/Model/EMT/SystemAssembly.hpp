/**
 * @file SystemAssembly.hpp
 *
 * @brief Free-function assembly utilities shared by EMT aggregates: contiguous
 * index layout, binding to the global arrays, and building or updating the
 * system CSR Jacobian from the element COO Jacobians.
 *
 */

#pragma once

#include <algorithm>
#include <numeric>
#include <vector>

#include <GridKit/LinearAlgebra/SparseMatrix/CsrMatrix.hpp>
#include <GridKit/MemoryUtilities/MemoryUtils.hpp>
#include <GridKit/Model/EMT/Element.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace SystemAssembly
    {
      template <typename ScalarT, typename IdxT>
      IdxT totalSize(const std::vector<Element<ScalarT, IdxT>*>& elements)
      {
        IdxT size = 0;
        for (const auto* element : elements)
        {
          size += element->size();
        }
        return size;
      }

      template <typename ScalarT, typename IdxT>
      void bindAll(const std::vector<Element<ScalarT, IdxT>*>& elements,
                   ScalarT*                                    y,
                   ScalarT*                                    yp,
                   ScalarT*                                    f,
                   IdxT*                                       index)
      {
        IdxT base = 0;
        for (auto* element : elements)
        {
          element->bind(y, yp, f, index, base);
          base += element->size();
        }
      }

      /**
       * @brief Build the system CSR Jacobian pattern from the element Jacobians
       * and a stable map from each gathered triplet to its CSR value slot.
       *
       * Element entries are already in global coordinates and each residual row
       * belongs to a single element, so the only duplicate (row, col) pairs come
       * from an element's own df/dy and alpha df/dy' contributions, which are
       * summed into one slot.
       */
      template <typename ScalarT, typename IdxT>
      void buildCsr(const std::vector<Element<ScalarT, IdxT>*>&                             elements,
                    IdxT                                                                    n,
                    LinearAlgebra::CsrMatrix<typename ScalarTraits<ScalarT>::RealT, IdxT>*& csr,
                    std::vector<IdxT>&                                                      map)
      {
        using RealT = typename ScalarTraits<ScalarT>::RealT;

        std::vector<IdxT> grow;
        std::vector<IdxT> gcol;
        for (auto* element : elements)
        {
          auto [rows, cols, vals] = element->jacobian().getEntries(false);
          for (size_t i = 0; i < rows.size(); ++i)
          {
            grow.push_back(rows[i]);
            gcol.push_back(cols[i]);
          }
        }

        const IdxT        total = static_cast<IdxT>(grow.size());
        std::vector<IdxT> order(static_cast<size_t>(total));
        std::iota(order.begin(), order.end(), IdxT{0});
        std::sort(order.begin(),
                  order.end(),
                  [&](IdxT a, IdxT b)
                  { return grow[a] < grow[b] || (grow[a] == grow[b] && gcol[a] < gcol[b]); });

        map.assign(static_cast<size_t>(total), 0);
        std::vector<IdxT> slot_row;
        std::vector<IdxT> slot_col;
        for (IdxT k = 0; k < total; ++k)
        {
          const IdxT idx = order[k];
          if (k == 0 || grow[idx] != slot_row.back() || gcol[idx] != slot_col.back())
          {
            slot_row.push_back(grow[idx]);
            slot_col.push_back(gcol[idx]);
          }
          map[idx] = static_cast<IdxT>(slot_col.size() - 1);
        }

        const IdxT nnz  = static_cast<IdxT>(slot_col.size());
        IdxT*      rows = new IdxT[static_cast<size_t>(n) + 1];
        IdxT*      cols = new IdxT[static_cast<size_t>(nnz)];
        RealT*     vals = new RealT[static_cast<size_t>(nnz)];

        for (IdxT i = 0; i <= n; ++i)
        {
          rows[i] = 0;
        }
        for (IdxT s = 0; s < nnz; ++s)
        {
          rows[slot_row[s] + 1]++;
          cols[s] = slot_col[s];
          vals[s] = 0.0;
        }
        for (IdxT i = 1; i <= n; ++i)
        {
          rows[i] += rows[i - 1];
        }

        csr = new LinearAlgebra::CsrMatrix<RealT, IdxT>(n, n, nnz, &rows, &cols, &vals);
      }

      template <typename ScalarT, typename IdxT>
      void updateCsrValues(const std::vector<Element<ScalarT, IdxT>*>&                            elements,
                           LinearAlgebra::CsrMatrix<typename ScalarTraits<ScalarT>::RealT, IdxT>* csr,
                           const std::vector<IdxT>&                                               map)
      {
        using RealT = typename ScalarTraits<ScalarT>::RealT;

        RealT* vals = csr->getValues();
        for (IdxT i = 0; i < csr->getNnz(); ++i)
        {
          vals[i] = 0.0;
        }

        IdxT k = 0;
        for (auto* element : elements)
        {
          auto [rows, cols, entries] = element->jacobian().getEntries(false);
          for (size_t i = 0; i < entries.size(); ++i)
          {
            vals[map[static_cast<size_t>(k)]] += entries[i];
            ++k;
          }
        }
        csr->setUpdated(GridKit::memory::HOST);
      }
    } // namespace SystemAssembly
  } // namespace EMT
} // namespace GridKit
