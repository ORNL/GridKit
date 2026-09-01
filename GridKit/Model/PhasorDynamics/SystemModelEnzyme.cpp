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
     * @brief Build the system Jacobian sparsity pattern and COO-to-CSR map.
     *
     * Runs once for a fixed topology, after every child has evaluated its
     * Jacobian and established its final COO structure.
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::buildJacobianStructure()
    {
      IdxT nnz_dup = 0;
      for (const auto* component : components_)
      {
        auto* component_jacobian = component->getCooJacobian();
        if (component_jacobian != nullptr)
        {
          nnz_dup += component_jacobian->getNnz();
        }
        else
        {
          Log::warning() << "A component has returned a nullptr Jacobian.\n";
        }
      }

      for (const auto* bus : buses_)
      {
        auto* bus_jacobian = bus->getCooJacobian();
        if (bus_jacobian != nullptr)
        {
          nnz_dup += bus_jacobian->getNnz();
        }
        else
        {
          Log::warning() << "A bus has returned a nullptr Jacobian.\n";
        }
      }

      IdxT*  rows_dup = new IdxT[static_cast<std::size_t>(nnz_dup)];
      IdxT*  cols_dup = new IdxT[static_cast<std::size_t>(nnz_dup)];
      RealT* vals_dup = new RealT[static_cast<std::size_t>(nnz_dup)];

      IdxT counter = 0;
      for (const auto* component : components_)
      {
        auto* component_jacobian = component->getCooJacobian();
        if (component_jacobian == nullptr)
        {
          continue;
        }

        const IdxT*  rows    = component_jacobian->getRowData();
        const IdxT*  columns = component_jacobian->getColData();
        const RealT* values  = component_jacobian->getValues();
        for (IdxT i = 0; i < component_jacobian->getNnz(); ++i, ++counter)
        {
          rows_dup[counter] = rows[i];
          cols_dup[counter] = columns[i];
          vals_dup[counter] = values[i];
        }
      }

      for (const auto* bus : buses_)
      {
        auto* bus_jacobian = bus->getCooJacobian();
        if (bus_jacobian == nullptr)
        {
          continue;
        }

        const IdxT*  rows    = bus_jacobian->getRowData();
        const IdxT*  columns = bus_jacobian->getColData();
        const RealT* values  = bus_jacobian->getValues();
        for (IdxT i = 0; i < bus_jacobian->getNnz(); ++i, ++counter)
        {
          rows_dup[counter] = rows[i];
          cols_dup[counter] = columns[i];
          vals_dup[counter] = values[i];
        }
      }

      CooMatrixT jac(size_, size_, nnz_dup, &rows_dup, &cols_dup, &vals_dup);
      IdxT*      row_ptrs = jac.getCsrRowData();

      nnz_ = jac.getNnz();

      IdxT*  cols = new IdxT[static_cast<std::size_t>(nnz_)];
      RealT* vals = new RealT[static_cast<std::size_t>(nnz_)];
      std::copy(jac.getColData(), jac.getColData() + nnz_, cols);
      std::copy(jac.getValues(), jac.getValues() + nnz_, vals);

      csr_jac_ = new CsrMatrixT(size_, size_, nnz_, &row_ptrs, &cols, &vals);

      const IdxT* map_to_sorted = jac.getMapToSorted();
      const IdxT* map_to_dedup  = jac.getMapToDeduplicated();
      map_to_csr_                = new IdxT[static_cast<std::size_t>(nnz_dup)];
      for (IdxT i = 0; i < nnz_dup; ++i)
      {
        map_to_csr_[map_to_sorted[i]] = map_to_dedup[i];
      }
    }

    /**
     * @brief Snapshot invariant Jacobian values and flatten varying entries.
     *
     * A component that supplies admittance stamps has a complete, constant
     * Jacobian contribution. Bus blocks are structural zeros. Everything else
     * is represented by stable pointers into its COO value buffer.
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::snapshotConstantJacobian()
    {
      constant_jacobian_values_.assign(static_cast<std::size_t>(nnz_), RealT{0});
      varying_jacobian_to_csr_.clear();
      varying_jacobian_sources_.clear();
      varying_jacobian_to_csr_.reserve(static_cast<std::size_t>(nnz_));
      varying_jacobian_sources_.reserve(static_cast<std::size_t>(nnz_));

      IdxT counter = 0;
      for (auto* component : components_)
      {
        auto* component_jacobian = component->getCooJacobian();
        if (component_jacobian == nullptr)
        {
          continue;
        }

        const bool   varies = component->admittanceStamps(nullptr) == 0;
        const RealT* values = component_jacobian->getValues();
        for (IdxT i = 0; i < component_jacobian->getNnz(); ++i, ++counter)
        {
          const IdxT destination = map_to_csr_[counter];
          if (varies)
          {
            varying_jacobian_to_csr_.push_back(destination);
            varying_jacobian_sources_.push_back(values + i);
          }
          else
          {
            constant_jacobian_values_[destination] += values[i];
          }
        }
      }

      // Finite-bus blocks reserve their diagonal structure with zero values.
      // Infinite buses have empty blocks, so neither kind contributes values.
      jacobian_snapshot_epoch_ = this->admittanceEpoch();
      jacobian_snapshot_ready_ = true;
    }

    /**
     * @brief Evaluate and assemble the system Jacobian.
     *
     * The CSR structure is built once per topology. Each call copies the
     * invariant stamped baseline and accumulates only entries from components
     * whose values may vary. An admittance epoch change forces a full child
     * evaluation and resnapshots both the values and stamped classification.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateJacobian()
    {
      const bool admittance_current = ensureAdmittanceCurrent();
      const bool refresh_snapshot   = csr_jac_ == nullptr
                                    || !jacobian_snapshot_ready_
                                    || jacobian_snapshot_epoch_ != this->admittanceEpoch()
                                    || !admittance_current;

      if (refresh_snapshot)
      {
        for (auto* bus : buses_)
        {
          bus->evaluateJacobian();
        }

        for (auto* component : components_)
        {
          component->evaluateJacobian();
        }

        if (csr_jac_ == nullptr)
        {
          buildJacobianStructure();
        }
        snapshotConstantJacobian();
      }
      else
      {
        for (auto* component : evaluated_components_)
        {
          component->evaluateJacobian();
        }
      }

      RealT* values = csr_jac_->getValues();
      std::copy(constant_jacobian_values_.begin(),
                constant_jacobian_values_.end(),
                values);

      for (std::size_t i = 0; i < varying_jacobian_to_csr_.size(); ++i)
      {
        values[varying_jacobian_to_csr_[i]] += *varying_jacobian_sources_[i];
      }

      return 0;
    }


    // Available template instantiations
    // template class SystemModel<double, long int>;
    template class SystemModel<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
