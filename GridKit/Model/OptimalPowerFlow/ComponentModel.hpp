/**
 * @file ComponentModel.hpp
 * @brief Evaluation of an optimal power flow component from its kernels.
 */

#pragma once

#include <string>
#include <vector>

#include <GridKit/Model/OptimalPowerFlow/Component.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /**
     * @brief Component evaluated from the kernels of `model_type`
     *
     * `model_type` provides `objective(x)` and `constraints(x, g)` over its
     * local variables, the sizes `VARIABLE_SIZE`, `CONSTRAINT_SIZE`,
     * `INTERNAL_SIZE` and `INTERNAL_CONSTRAINT_SIZE`, and `CONSTANT`, which is
     * true when neither kernel depends on the variables. `g` starts at zero,
     * so a kernel writes only the rows it contributes to. Derivatives come
     * from Enzyme for double scalars and from the system dependencies for
     * `DependencyTracking::Variable`.
     *
     * @tparam model_type - Derived component model
     */
    template <typename model_type, typename scalar_type, typename index_type>
    class ComponentModel : public Component<scalar_type, index_type>
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename Component<ScalarT, IdxT>::RealT;

      int evaluateObjective(const ScalarT* x, ScalarT& f) override;
      int evaluateGradient(const ScalarT* x, RealT* gradient) override;
      int evaluateConstraints(const ScalarT* x, ScalarT* g) override;
      int evaluateJacobian(const ScalarT* x) override;
      int evaluateHessian(const ScalarT* x, RealT sigma, const RealT* lambda) override;
      int evaluateTerminalPower(const ScalarT* x, ScalarT* power) override;

    protected:
      /**
       * @param[in] id - Case device `id`
       * @param[in] terminals - Bus number of each terminal
       */
      ComponentModel(std::string id, std::vector<IdxT> terminals);

    private:
      static constexpr bool hasEnzymeDerivatives();

      const model_type& model() const;

      auto gather(const ScalarT* x) const;
    };
  } // namespace OptimalPowerFlow
} // namespace GridKit
