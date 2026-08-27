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
     *
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateJacobian()
    {
      using ProfileClock = std::chrono::steady_clock;

      const auto jacobian_start = ProfileClock::now();

      // Initialize bus Jacobians
      const auto bus_start = ProfileClock::now();
      for (const auto& bus : buses_)
      {
        bus->evaluateJacobian();
      }
      profile_bus_jacobian_seconds_ += std::chrono::duration<double>(ProfileClock::now() - bus_start).count();

      // Evaluate component Jacobians, including contribution to the bus Jacobians
      IdxT profile_begin = 0;
      for (std::size_t group = 0; group < profile_group_count_; ++group)
      {
        const IdxT profile_end = profile_component_ends_[group];
        if (profile_begin < profile_end)
        {
          const auto group_start = ProfileClock::now();
          for (IdxT i = profile_begin; i < profile_end; ++i)
          {
            components_[static_cast<std::size_t>(i)]->evaluateJacobian();
          }
          profile_jacobian_seconds_[group] += std::chrono::duration<double>(ProfileClock::now() - group_start).count();
        }
        profile_begin = profile_end;
      }

      const IdxT component_count = static_cast<IdxT>(components_.size());
      if (profile_begin < component_count)
      {
        const auto other_start = ProfileClock::now();
        for (IdxT i = profile_begin; i < component_count; ++i)
        {
          components_[static_cast<std::size_t>(i)]->evaluateJacobian();
        }
        profile_other_jacobian_seconds_ += std::chrono::duration<double>(ProfileClock::now() - other_start).count();
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

        // Build the system COO Jacobian
        CooMatrixT jac(size_, size_, nnz_dup, &rows_dup, &cols_dup, &vals_dup);

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

        // Build a mappping from original COO index to CSR index
        map_to_csr_ = new IdxT[static_cast<size_t>(nnz_dup)];
        for (IdxT i = 0; i < nnz_dup; ++i)
        {
          map_to_csr_[map_to_sorted[i]] = map_to_dedup[i];
        }
      }
      else
      {
        // Zero out values
        const auto zero_start = ProfileClock::now();
        RealT*     vals       = csr_jac_->getValues();
        for (IdxT i = 0; i < csr_jac_->getNnz(); ++i)
        {
          vals[i] = 0.0;
        }
        profile_jacobian_zero_seconds_ += std::chrono::duration<double>(ProfileClock::now() - zero_start).count();

        // Update CSR values from component and bus Jacobians
        const auto component_scatter_start = ProfileClock::now();
        IdxT       counter                 = 0;
        for (const auto& component : components_)
        {
          auto component_jacobian = component->getCooJacobian();

          if (component_jacobian != nullptr)
          {
            const RealT* values = component_jacobian->getValues();
            for (IdxT i = 0; i < component_jacobian->getNnz(); ++i)
            {
              vals[map_to_csr_[counter]] += values[i];
              counter++;
            }
          }
        }
        profile_component_scatter_seconds_ += std::chrono::duration<double>(ProfileClock::now() - component_scatter_start).count();

        const auto bus_scatter_start = ProfileClock::now();
        for (const auto& bus : buses_)
        {
          auto bus_jacobian = bus->getCooJacobian();

          if (bus_jacobian != nullptr)
          {
            const RealT* values = bus_jacobian->getValues();
            for (IdxT i = 0; i < bus_jacobian->getNnz(); ++i)
            {
              vals[map_to_csr_[counter]] += values[i];
              counter++;
            }
          }
        }
        profile_bus_scatter_seconds_ += std::chrono::duration<double>(ProfileClock::now() - bus_scatter_start).count();
      }

      // Log::misc() << "System Enzyme Jacobian\n";
      // csr_jac_->print(Log::misc());

      profile_system_jacobian_seconds_ += std::chrono::duration<double>(ProfileClock::now() - jacobian_start).count();
      ++profile_jacobian_calls_;

      return 0;
    }

    // Available template instantiations
    // template class SystemModel<double, long int>;
    template class SystemModel<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
