/**
 * @brief A collection of mixins meant to provide common functionality in implementing the
 * `GridKit::Model::Evaluator` interface.
 */
#pragma once

#include <optional>

#include "Evaluator.hpp"

namespace GridKit
{
  namespace Mixin
  {
    namespace Evaluator
    {
      template <class ScalarT, typename IdxT, class Derived>
      class CsrJacobian : public virtual Model::Evaluator<ScalarT, IdxT>
      {
        using RealT = typename Model::Evaluator<ScalarT, IdxT>::RealT;
        using JacT  = typename Model::Evaluator<ScalarT, IdxT>::CsrJacobian;

      protected:
        /**
         * @brief The jacobian matrix.
         * @note The jacobian isn't computed (and doesn't exist) until `evaluateJacobian()` is called.
         */
        JacT csr_jacobian_ = JacT(Derived::SIZE, Derived::SIZE);

      public:
        bool hasCsrJacobian() final
        {
          return true;
        }

        int evaluateCsrJacobian() override
        {
          using CsrBuilder = LinearAlgebra::CsrBuilder<RealT, IdxT>;

          auto& self = static_cast<Derived&>(*this);

          if (csr_jacobian_.numNonZero() == 0)
          {
            csr_jacobian_ = self.buildCsrJacobian(CsrBuilder::fromEmpty(Derived::SIZE, Derived::SIZE));
          }
          else
          {
            csr_jacobian_ = self.buildCsrJacobian(CsrBuilder::fromTemplate(std::move(csr_jacobian_)));
          }

          return 0;
        }

        JacT& getCsrJacobian() final
        {
          return csr_jacobian_;
        }

        const JacT& getCsrJacobian() const final
        {
          return csr_jacobian_;
        }
      };
    } // namespace Evaluator
  } // namespace Mixin
} // namespace GridKit