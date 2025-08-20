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
      class CSRJacobian : public virtual Model::Evaluator<ScalarT, IdxT>
      {
        using JacT = typename Model::Evaluator<ScalarT, IdxT>::CSRJacobian;

      protected:
        /**
         * @brief The jacobian matrix.
         * @note The jacobian isn't computed (and doesn't exist) until `evaluateJacobian()` is called.
         */
        JacT csr_jacobian_ = JacT(Derived::SIZE, Derived::SIZE);

      public:
        bool hasCSRJacobian() final
        {
          return true;
        }

        int evaluateCSRJacobian() override
        {
          using CSRBuilder = CSRBuilder<ScalarT, IdxT>;

          auto build = [](auto builder)
          {
            builder.row(0).elem(0, 1.0 / R_).elem(1, -1.0 / R_);
            builder.row(1).elem(0, -1.0 / R_).elem(1, 1.0 / R_);

            return builder;
          };

          if (csr_jacobian_.numNonZero() == 0)
          {
            csr_jacobian_ = Derived::buildCSR(CSRBuilder::fromEmpty(Derived::SIZE, Derived::SIZE));
          }
          else
          {
            csr_jacobian_ = Derived::buildCSR(CSRBuilder::fromTemplate(std::move(csr_jacobian_)));
          }
        }

        JacT& getCSRJacobian() final
        {
          return csr_jacobian_;
        }

        const JacT& getCSRJacobian() final const
        {
          return csr_jacobian_;
        }
      };
    } // namespace Evaluator
  } // namespace Mixin
} // namespace GridKit