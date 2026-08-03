/**
 * @file ShuntReduction.cpp
 *
 * @brief Initialization, equation residual, and sparse-jet Jacobian for the
 * shunt reduction.
 *
 */

#include "ShuntReduction.hpp"

#include <utility>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename ScalarT, typename IdxT>
      ShuntReduction<ScalarT, IdxT>::ShuntReduction(
          const ShuntView<ScalarT, IdxT>& full, MapT map)
        : full_(full),
          map_(std::move(map)),
          layout_(map_.conductors, map_.terminals)
      {
      }

      template <typename ScalarT, typename IdxT>
      int ShuntReduction<ScalarT, IdxT>::initialize()
      {
        const Layout& layout = layout_;
        auto          G      = input(full_.G());
        auto          C      = input(full_.C());
        auto          Gred   = state(layout.g, layout.P, layout.P);
        auto          Cred   = state(layout.c, layout.P, layout.P);

        for (IdxT t = 0; t < layout.P; ++t)
        {
          for (IdxT q = 0; q < layout.P; ++q)
          {
            ScalarT sum_g{};
            ScalarT sum_c{};
            for (const IdxT a : map_.members[static_cast<size_t>(t)])
            {
              for (const IdxT b : map_.members[static_cast<size_t>(q)])
              {
                sum_g += G(a, b);
                sum_c += C(a, b);
              }
            }
            Gred(t, q) = sum_g;
            Cred(t, q) = sum_c;
          }
        }
        return 0;
      }

      template <typename ScalarT, typename IdxT>
      template <typename V>
      __attribute__((always_inline)) inline void
      ShuntReduction<ScalarT, IdxT>::evaluateResidualKernel(V* y,
                                                            V*,
                                                            V* u,
                                                            V* f) const
      {
        const Layout& layout = layout_;
        auto          Gred   = slice(y, layout.g, layout.P, layout.P);
        auto          Cred   = slice(y, layout.c, layout.P, layout.P);
        auto          FG     = slice(f, layout.g, layout.P, layout.P);
        auto          FC     = slice(f, layout.c, layout.P, layout.P);
        auto          G      = input(u, full_.G());
        auto          C      = input(u, full_.C());

        // -Yred + E^T Y E = 0: pure gathers through the incidence.
        for (IdxT t = 0; t < layout.P; ++t)
        {
          for (IdxT q = 0; q < layout.P; ++q)
          {
            V sum_g{};
            V sum_c{};
            for (const IdxT a : map_.members[static_cast<size_t>(t)])
            {
              for (const IdxT b : map_.members[static_cast<size_t>(q)])
              {
                sum_g += G(a, b);
                sum_c += C(a, b);
              }
            }
            FG(t, q) = sum_g - Gred(t, q);
            FC(t, q) = sum_c - Cred(t, q);
          }
        }
      }

      template <typename ScalarT, typename IdxT>
      void ShuntReduction<ScalarT, IdxT>::evaluateInternalResidual(
          ScalarT* y, ScalarT* yp, ScalarT* u, ScalarT* f) const
      {
        evaluateResidualKernel<ScalarT>(y, yp, u, f);
      }

      template <typename ScalarT, typename IdxT>
      int ShuntReduction<ScalarT, IdxT>::evaluateResidual()
      {
        evaluateInternalResidual(y_, yp_, nullptr, f_);
        return 0;
      }

      template <typename ScalarT, typename IdxT>
      int ShuntReduction<ScalarT, IdxT>::evaluateJacobian()
      {
        return this->template evaluateInputAlgebraicJacobian<
            ShuntReduction<ScalarT, IdxT>>();
      }

      template class ShuntReduction<double, long int>;
      template class ShuntReduction<double, size_t>;
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
