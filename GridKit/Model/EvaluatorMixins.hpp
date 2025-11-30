/**
 * @brief A collection of mixins meant to provide common functionality in implementing the
 * `GridKit::Model::Evaluator` interface.
 */
#pragma once

#include <concepts>
#include <optional>
#include <utility>

#include "Evaluator.hpp"

namespace GridKit
{
  namespace Mixin
  {
    namespace Evaluator
    {
      /**
       * @brief Describes types which may use `GridKit::Mixin::Evaluator::CsrJacobian` as a mixin.
       */
      template <class Derived, class CsrBuilder, class CsrJacobian>
      concept CanMixinCsrJacobian = requires(Derived d, CsrJacobian j) {
        { Derived::SIZE } -> std::convertible_to<size_t>;
        { d.buildCsrJacobian(CsrBuilder::fromEmpty(Derived::SIZE, Derived::SIZE)) } -> std::convertible_to<CsrJacobian>;
        { d.buildCsrJacobian(CsrBuilder::fromTemplate(std::move(j))) } -> std::convertible_to<CsrJacobian>;
      };

      /**
       * @brief  Implements the common pattern of how many components will want to compute their CSR Jacobian.
       *         Store a square CSR matrix based on the component size, and the first time `evaluateJacobian()`
       *         is called, use `LinearAlgebra::CsrBuilder::fromEmpty()`. After that, use `LinearAlgebra::CsrBuilder::fromTemplate`, using the
       *         previously computed jacobian as a template.
       * @tparam Derived The class which is inheriting this mixin. Must meet the equirements of `CanMixinCsrJacobian`.
       *
       * An example of how a component might use the mixin:
       * @code{.cpp}
       * template <class ScalarT, typename IdxT>
       * class Component : public virtual Evaluator<ScalarT, IdxT>, public Mixin::Evaluator::CsrJacobian<ScalarT, IdxT, Resistor<ScalarT, IdxT>>
       * {
       *  public:
       *    const static size_t SIZE = 2;
       *
       *    template <bool INCLUDE_DIAGONALS, bool KEEP_SORTED, bool USE_TEMPLATE>
       *    CsrJacobian buildCsrJacobian(LinearAlgebra::CsrBuilder<ScalarT, IdxT, INCLUDE_DIAGONALS, KEEP_SORTED, USE_TEMPLATE> builder);
       * };
       * @endcode
       *
       * @note   Assumes that if `::csr_jacobian_` is non-empty that the next call to `buildCsrJacobian()` will have the same
       *         sparsity pattern as the last call. If this is not true, this mixin cannot be used correctly.
       * @note   Since `CsrJacobian` only partially implements `Model::Evaluator`, any component which wishes to use this mixin
       *         must inherit `Model::Evaluator` virtually, and implement the remaining members.
       */
      template <class ScalarT, typename IdxT, template <class S, typename I> class Derived>
      class CsrJacobian : public virtual Model::Evaluator<ScalarT, IdxT>
      {
        using RealT    = typename Model::Evaluator<ScalarT, IdxT>::RealT;
        using JacT     = typename Model::Evaluator<ScalarT, IdxT>::CsrJacobian;
        using DerivedT = Derived<ScalarT, IdxT>;

      protected:
        /**
         * @brief The jacobian matrix.
         * @note The jacobian isn't computed (and is empty) until `evaluateJacobian()` is called.
         */
        JacT csr_jacobian_ = JacT(DerivedT::SIZE, DerivedT::SIZE);

        bool should_use_csr_ = true;

        CsrJacobian(const CsrJacobian& other) : should_use_csr_(other.should_use_csr_)
        {
        }

        CsrJacobian() = default;

      public:
        bool hasCsrJacobian() final
        {
          return should_use_csr_;
        }

        /**
         * @brief Evaluate the jacobian using a common pattern for components - if the jacobian is empty, use
         * `LinearAlgebra::CsrBuilder::fromEmpty()` and if the jacobian has been previously computed, use it as
         * a template. Calls `Derived::buildCsrJacobian()` with a builder, and expects something convertible to
         * `JacT` as a returned value for the value of the Jacobian.
         */
        int evaluateCsrJacobian() override
        {
          using CsrBuilder = LinearAlgebra::CsrBuilder<RealT, IdxT>;

          static_assert(CanMixinCsrJacobian<DerivedT, CsrBuilder, JacT>);

          auto& self = static_cast<DerivedT&>(*this);

          if (csr_jacobian_.numNonZero() == 0)
          {
            csr_jacobian_ = self.buildCsrJacobian(CsrBuilder::fromEmpty(DerivedT::SIZE, DerivedT::SIZE));
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

        void setCsrUsage(bool use) noexcept
        {
          should_use_csr_ = use;
        }
      };
    } // namespace Evaluator
  } // namespace Mixin
} // namespace GridKit
