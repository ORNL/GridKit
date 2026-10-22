#include "SystemModelImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Check components for Jacobian availability
     *
     * @return true
     * @return false
     */
    template <typename scalar_type, typename index_type>
    bool SystemModel<scalar_type, index_type>::hasJacobian()
    {
      bool has_jacobian = true;
      for (const auto& component : components_)
      {
        has_jacobian = has_jacobian && component->hasJacobian();
      }

      for (const auto& bus : buses_)
      {
        has_jacobian = has_jacobian && bus->hasJacobian();
      }

      if (!has_jacobian)
      {
        Log::warning() << "GridKit was built with Enzyme, but some models "
                          "don't have Jacobians available. "
                          "Falling back to dense Jacobians in PhasorDynamics.\n";
      }

      return has_jacobian;
    }

    /**
     * @brief Evaluate system Jacobian using component-level Jacobians.
     *
     * - Evaluate component-level Jacobians (stored as CooMatrixT).
     * - If not already constructed, construct system-level Jacobian.
     *   This used a system-level CooMatrix to deduplicate and sort the sparsity pattern.
     * - If already constructed, reset Jacobian values to zero and accumulate contributions.
     * - Entries of computed outputs enter through the chain rule,
     *   \f$ J_{xx} + J_{xo} J_{ox} \f$, never as solver variables.
     *
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateJacobian()
    {
      // Initialize bus Jacobians
      for (const auto& bus : buses_)
      {
        bus->evaluateJacobian();
      }

      // Evaluate component Jacobians, including contribution to the bus Jacobians
      for (const auto& component : components_)
      {
        component->evaluateJacobian();
      }

      // Build or update system CSR Jacobian
      if (csr_jac_ == nullptr)
      {
        // Count the number of non-zeros
        IdxT nnz_dup = 0;
        for (const auto& component : components_)
        {
          const auto component_jacobian = component->getCooJacobian();

          if (component_jacobian != nullptr)
          {
            nnz_dup += component_jacobian->getNnz();
          }
          else
          {
            Log::warning() << "A component has returned a nullptr Jacobian.\n";
          }
        }

        for (const auto& bus : buses_)
        {
          auto bus_jacobian = bus->getCooJacobian();

          if (bus_jacobian != nullptr)
          {
            nnz_dup += bus_jacobian->getNnz();
          }
          else
          {
            Log::warning() << "A bus has returned a nullptr Jacobian.\n";
          }
        }

        // Allocate COO triplet arrays (we own these until we hand off to CsrMatrix)
        IdxT*  rows_dup = new IdxT[static_cast<size_t>(nnz_dup)];
        IdxT*  cols_dup = new IdxT[static_cast<size_t>(nnz_dup)];
        RealT* vals_dup = new RealT[static_cast<size_t>(nnz_dup)];

        IdxT counter = 0;
        for (const auto& component : components_)
        {
          auto component_jacobian = component->getCooJacobian();

          if (component_jacobian != nullptr)
          {
            const IdxT*  rows    = component_jacobian->getRowData();
            const IdxT*  columns = component_jacobian->getColData();
            const RealT* values  = component_jacobian->getValues();
            for (IdxT i = 0; i < component_jacobian->getNnz(); ++i)
            {
              rows_dup[counter] = rows[i];
              cols_dup[counter] = columns[i];
              vals_dup[counter] = values[i];
              counter++;
            }
          }
          else
          {
            Log::warning() << "A component has returned a nullptr Jacobian.\n";
          }
        }

        for (const auto& bus : buses_)
        {
          auto bus_jacobian = bus->getCooJacobian();

          if (bus_jacobian != nullptr)
          {
            const IdxT*  rows    = bus_jacobian->getRowData();
            const IdxT*  columns = bus_jacobian->getColData();
            const RealT* values  = bus_jacobian->getValues();
            for (IdxT i = 0; i < bus_jacobian->getNnz(); ++i)
            {
              rows_dup[counter] = rows[i];
              cols_dup[counter] = columns[i];
              vals_dup[counter] = values[i];
              counter++;
            }
          }
          else
          {
            Log::warning() << "A bus has returned a nullptr Jacobian.\n";
          }
        }

        // Computed outputs are not variables. Entries in their rows are output
        // gradients and entries in their columns are consumer partials; the
        // chain rule replaces each consumer partial with its products through
        // the gradient of the output it reads.
        std::vector<std::vector<IdxT>> output_gradients(static_cast<size_t>(output_count_));
        IdxT                           nnz_direct = 0;
        for (IdxT i = 0; i < nnz_dup; ++i)
        {
          const bool output_row    = rows_dup[i] >= size_;
          const bool output_column = cols_dup[i] >= size_;

          if (output_row && output_column)
          {
            Log::error() << "SystemModel: a computed output cannot depend on another computed output\n";
            delete[] rows_dup;
            delete[] cols_dup;
            delete[] vals_dup;
            return 1;
          }

          if (output_row)
          {
            output_gradients[static_cast<size_t>(rows_dup[i] - size_)].push_back(i);
          }

          if (!output_row && !output_column)
          {
            ++nnz_direct;
          }
        }

        chain_entries_.clear();
        for (IdxT i = 0; i < nnz_dup; ++i)
        {
          if (rows_dup[i] < size_ && cols_dup[i] >= size_)
          {
            for (const IdxT j : output_gradients[static_cast<size_t>(cols_dup[i] - size_)])
            {
              chain_entries_.push_back({i, j, INVALID_INDEX<IdxT>});
            }
          }
        }

        // Direct entries followed by chain-rule products
        const IdxT nnz_coo  = nnz_direct + static_cast<IdxT>(chain_entries_.size());
        IdxT*      rows_coo = new IdxT[static_cast<size_t>(nnz_coo)];
        IdxT*      cols_coo = new IdxT[static_cast<size_t>(nnz_coo)];
        RealT*     vals_coo = new RealT[static_cast<size_t>(nnz_coo)];

        std::vector<IdxT> coo_entry(static_cast<size_t>(nnz_dup), INVALID_INDEX<IdxT>);
        IdxT              k = 0;
        for (IdxT i = 0; i < nnz_dup; ++i)
        {
          if (rows_dup[i] < size_ && cols_dup[i] < size_)
          {
            rows_coo[k]  = rows_dup[i];
            cols_coo[k]  = cols_dup[i];
            vals_coo[k]  = vals_dup[i];
            coo_entry[i] = k;
            ++k;
          }
        }

        for (const auto& entry : chain_entries_)
        {
          rows_coo[k] = rows_dup[entry.consumer];
          cols_coo[k] = cols_dup[entry.output];
          vals_coo[k] = vals_dup[entry.consumer] * vals_dup[entry.output];
          ++k;
        }

        output_entry_values_.clear();
        if (nnz_direct < nnz_dup)
        {
          output_entry_values_.assign(static_cast<size_t>(nnz_dup), RealT{0.0});
        }

        delete[] rows_dup;
        delete[] cols_dup;
        delete[] vals_dup;

        // Build the system COO Jacobian
        CooMatrixT jac(size_, size_, nnz_coo, &rows_coo, &cols_coo, &vals_coo);

        // Populate CSR data with sort and deduplicate
        IdxT* row_ptrs = jac.getCsrRowData();

        // Deduplicated nnz
        nnz_ = jac.getNnz();

        // Allocate cols/vals with deduplicated nnz
        IdxT*  cols = new IdxT[static_cast<size_t>(nnz_)];
        RealT* vals = new RealT[static_cast<size_t>(nnz_)];

        std::copy(jac.getColData(), jac.getColData() + nnz_, cols);
        std::copy(jac.getValues(), jac.getValues() + nnz_, vals);

        // Create the CSR Jacobian
        csr_jac_ = new CsrMatrixT(size_, size_, nnz_, &row_ptrs, &cols, &vals);

        const IdxT* map_to_sorted = jac.getMapToSorted();
        const IdxT* map_to_dedup  = jac.getMapToDeduplicated();

        // CSR slot of each system COO entry
        std::vector<IdxT> coo_to_csr(static_cast<size_t>(nnz_coo));
        for (IdxT i = 0; i < nnz_coo; ++i)
        {
          coo_to_csr[map_to_sorted[i]] = map_to_dedup[i];
        }

        // Build a mapping from original COO index to CSR index; entries in
        // computed-output rows or columns reach the CSR through chain_entries_
        map_to_csr_ = new IdxT[static_cast<size_t>(nnz_dup)];
        for (IdxT i = 0; i < nnz_dup; ++i)
        {
          map_to_csr_[i] = INVALID_INDEX<IdxT>;
          if (coo_entry[i] != INVALID_INDEX<IdxT>)
          {
            map_to_csr_[i] = coo_to_csr[coo_entry[i]];
          }
        }

        for (size_t c = 0; c < chain_entries_.size(); ++c)
        {
          chain_entries_[c].slot = coo_to_csr[static_cast<size_t>(nnz_direct) + c];
        }
      }
      else
      {
        // Zero out values
        RealT* vals = csr_jac_->getValues();
        for (IdxT i = 0; i < csr_jac_->getNnz(); ++i)
        {
          vals[i] = 0.0;
        }

        // Update CSR values from component and bus Jacobians
        IdxT counter = 0;
        for (const auto& component : components_)
        {
          auto component_jacobian = component->getCooJacobian();

          if (component_jacobian != nullptr)
          {
            const RealT* values = component_jacobian->getValues();
            for (IdxT i = 0; i < component_jacobian->getNnz(); ++i)
            {
              const IdxT slot = map_to_csr_[counter];
              if (slot != INVALID_INDEX<IdxT>)
              {
                vals[slot] += values[i];
              }
              else
              {
                output_entry_values_[static_cast<size_t>(counter)] = values[i];
              }
              counter++;
            }
          }
        }

        for (const auto& bus : buses_)
        {
          auto bus_jacobian = bus->getCooJacobian();

          if (bus_jacobian != nullptr)
          {
            const RealT* values = bus_jacobian->getValues();
            for (IdxT i = 0; i < bus_jacobian->getNnz(); ++i)
            {
              const IdxT slot = map_to_csr_[counter];
              if (slot != INVALID_INDEX<IdxT>)
              {
                vals[slot] += values[i];
              }
              else
              {
                output_entry_values_[static_cast<size_t>(counter)] = values[i];
              }
              counter++;
            }
          }
        }

        // Chain rule for computed outputs
        for (const auto& entry : chain_entries_)
        {
          vals[entry.slot] += output_entry_values_[static_cast<size_t>(entry.consumer)]
                              * output_entry_values_[static_cast<size_t>(entry.output)];
        }
      }

      // Log::misc() << "System Enzyme Jacobian\n";
      // csr_jac_->print(Log::misc());

      return 0;
    }

    // Available template instantiations
    // template class SystemModel<double, long int>;
    template class SystemModel<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
