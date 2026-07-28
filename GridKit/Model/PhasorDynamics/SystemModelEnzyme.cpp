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
      if (csr_jac_ == nullptr)
      {
        for (const auto& bus : buses_)
        {
          bus->evaluateJacobian();
        }

        for (const auto& component : components_)
        {
          component->evaluateJacobian();
        }

        buildJacobianStructure();
        snapshotConstantJacobian();
      }

      for (const auto& component : evaluated_components_)
      {
        component->evaluateJacobian();
      }

      RealT* vals = csr_jac_->getValues();
      std::copy(constant_values_.begin(), constant_values_.end(), vals);

      const std::size_t entries = block_to_csr_.size();
      for (std::size_t i = 0; i < entries; ++i)
      {
        vals[block_to_csr_[i]] += *block_source_[i];
      }

      return 0;
    }

    // Available template instantiations
    // template class SystemModel<double, long int>;
    template class SystemModel<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
