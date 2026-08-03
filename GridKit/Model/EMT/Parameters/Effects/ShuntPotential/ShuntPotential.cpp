/**
 * @file ShuntPotential.cpp
 *
 * @brief Initialization, equation residual, and sparse-jet Jacobian for the ShuntPotential model.
 *
 */

#include "ShuntPotential.hpp"

#include <algorithm>
#include <cassert>
#include <limits>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename ScalarT, typename IdxT>
      ShuntPotential<ScalarT, IdxT>::ShuntPotential(const Tower<ScalarT, IdxT>&     tower,
                                                    const Conductor<ScalarT, IdxT>& conductor)
        : tower_(tower),
          conductor_(conductor),
          layout_(conductor.conductorCount())
      {
      }

      template <typename ScalarT, typename IdxT>
      template <typename V>
      __attribute__((always_inline)) inline void
      ShuntPotential<ScalarT, IdxT>::evaluateResidualKernel(V* y,
                                                            V*,
                                                            V*,
                                                            V* f) const
      {
        using namespace GridKit::LinearAlgebra;

        const Layout& L      = layout_;
        const auto    scale  = 2.0 * Constants::pi<ScalarT>() * Constants::epsilon0<ScalarT>();
        auto          P      = slice(y, L.p, L.K, L.K);
        auto          G      = slice(y, L.g, L.K, L.K);
        auto          C      = slice(y, L.c, L.K, L.K);
        auto          Pres   = slice(f, L.p, L.K, L.K);
        auto          Gres   = slice(f, L.g, L.K, L.K);
        auto          Cres   = slice(f, L.c, L.K, L.K);
        auto          Lambda = matrix(L.K,
                             L.K,
                             [&](IdxT i, IdxT j)
                             {
                               return lambda(i, j);
                             });

        equation(Pres) = -scale * P + Lambda;
        equation(Gres) = -G;
        equation(Cres) = P * C - identity(P);
      }

      template <typename ScalarT, typename IdxT>
      void ShuntPotential<ScalarT, IdxT>::evaluateInternalResidual(ScalarT* y,
                                                                   ScalarT* yp,
                                                                   ScalarT* u,
                                                                   ScalarT* f) const
      {
        evaluateResidualKernel<ScalarT>(y, yp, u, f);
      }

      template <typename ScalarT, typename IdxT>
      int ShuntPotential<ScalarT, IdxT>::evaluateResidual()
      {
        evaluateInternalResidual(y_, yp_, nullptr, f_);
        return 0;
      }

      template <typename ScalarT, typename IdxT>
      int ShuntPotential<ScalarT, IdxT>::initialize()
      {
        const Layout& L     = layout_;
        const ScalarT scale = 2.0 * Constants::pi<ScalarT>() * Constants::epsilon0<ScalarT>();
        auto          P     = state(L.p, L.K, L.K);
        auto          G     = state(L.g, L.K, L.K);

        for (IdxT i = 0; i < L.K; ++i)
        {
          for (IdxT j = 0; j < L.K; ++j)
          {
            P(i, j) = lambda(i, j) / scale;
            G(i, j) = 0.0;
          }
        }

        solveCapacitance(y_);
        return 0;
      }

      template <typename ScalarT, typename IdxT>
      void ShuntPotential<ScalarT, IdxT>::solveCapacitance(ScalarT* y) const
      {
        const Layout&        L = layout_;
        std::vector<ScalarT> a(static_cast<size_t>(L.KK), 0.0);
        std::vector<ScalarT> b(static_cast<size_t>(L.KK), 0.0);
        Slice<ScalarT, IdxT> A{a.data(), nullptr, L.K, L.K};
        Slice<ScalarT, IdxT> B{b.data(), nullptr, L.K, L.K};
        auto                 P = slice(y, L.p, L.K, L.K);
        auto                 C = slice(y, L.c, L.K, L.K);

        for (IdxT i = 0; i < L.K; ++i)
        {
          for (IdxT j = 0; j < L.K; ++j)
          {
            A(i, j) = P(i, j);
          }
          B(i, i) = 1.0;
        }

        for (IdxT k = 0; k < L.K; ++k)
        {
          IdxT    pivot    = k;
          ScalarT max_norm = std::abs(A(k, k));
          for (IdxT i = k + 1; i < L.K; ++i)
          {
            const ScalarT row_norm = std::abs(A(i, k));
            if (row_norm > max_norm)
            {
              pivot    = i;
              max_norm = row_norm;
            }
          }

          assert(max_norm > std::numeric_limits<ScalarT>::epsilon());

          if (pivot != k)
          {
            for (IdxT j = 0; j < L.K; ++j)
            {
              std::swap(A(k, j), A(pivot, j));
              std::swap(B(k, j), B(pivot, j));
            }
          }

          const ScalarT diag = A(k, k);
          for (IdxT i = k + 1; i < L.K; ++i)
          {
            const ScalarT factor = A(i, k) / diag;
            A(i, k)              = 0.0;
            for (IdxT j = k + 1; j < L.K; ++j)
            {
              A(i, j) -= factor * A(k, j);
            }
            for (IdxT q = 0; q < L.K; ++q)
            {
              B(i, q) -= factor * B(k, q);
            }
          }
        }

        for (IdxT q = 0; q < L.K; ++q)
        {
          for (IdxT ii = L.K; ii > 0; --ii)
          {
            const IdxT i   = ii - 1;
            ScalarT    sum = B(i, q);
            for (IdxT j = i + 1; j < L.K; ++j)
            {
              sum -= A(i, j) * B(j, q);
            }
            B(i, q) = sum / A(i, i);
          }
        }

        for (IdxT i = 0; i < L.K; ++i)
        {
          for (IdxT q = 0; q < L.K; ++q)
          {
            C(i, q) = B(i, q);
          }
        }
      }

      template <typename ScalarT, typename IdxT>
      int ShuntPotential<ScalarT, IdxT>::evaluateJacobian()
      {
        return this->template evaluateLocalAlgebraicJacobian<ShuntPotential<ScalarT, IdxT>>();
      }

      template class ShuntPotential<double, long int>;
      template class ShuntPotential<double, size_t>;
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
