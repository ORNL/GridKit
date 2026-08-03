/**
 * @file SkinEffect.cpp
 *
 * @brief Initialization, equation residual, and sparse-jet Jacobian for the SkinEffect model.
 *
 */

#include "SkinEffect.hpp"

#include <GridKit/LinearAlgebra/DenseExpression.hpp>
#include <GridKit/Model/EMT/Constants.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename ScalarT, typename IdxT>
      SkinEffect<ScalarT, IdxT>::SkinEffect(const Conductor<ScalarT, IdxT>& conductor)
        : conductor_(conductor),
          layout_(conductor.conductorCount(), rootCount())
      {
      }

      template <typename ScalarT, typename IdxT>
      ScalarT SkinEffect<ScalarT, IdxT>::branchResistance(IdxT i, IdxT k) const
      {
        const ScalarT pi    = Constants::pi<ScalarT>();
        const ScalarT r     = conductor_.radius(i);
        const ScalarT sigma = conductor_.sigma(i);
        const ScalarT xi    = (static_cast<ScalarT>(k) + 0.75) * pi;
        return xi * xi / (4.0 * pi * sigma * r * r);
      }

      template <typename ScalarT, typename IdxT>
      ScalarT SkinEffect<ScalarT, IdxT>::branchInductance(IdxT i) const
      {
        return conductor_.mu(i) / (4.0 * Constants::pi<ScalarT>());
      }

      template <typename ScalarT, typename IdxT>
      int SkinEffect<ScalarT, IdxT>::initialize()
      {
        const Layout& L     = layout_;
        const ScalarT omega = static_cast<ScalarT>(coordinate());
        auto          Q     = state(L.q, L.K, L.S);
        auto          G     = state(L.g, L.K);
        auto          H     = state(L.h, L.K);
        auto          W     = state(L.w, L.K);
        auto          R     = state(L.r, L.K);
        auto          Li    = state(L.l, L.K);

        for (IdxT i = 0; i < L.K; ++i)
        {
          const ScalarT lb = branchInductance(i);
          G[i]             = 0.0;
          H[i]             = 0.0;
          for (IdxT k = 0; k < L.S; ++k)
          {
            const ScalarT rb  = branchResistance(i, k);
            Q(i, k)           = rb * rb + omega * omega * lb * lb;
            G[i]             += rb / Q(i, k);
            H[i]             += lb / Q(i, k);
          }
          W[i]  = G[i] * G[i] + omega * omega * H[i] * H[i];
          R[i]  = G[i] / W[i];
          Li[i] = H[i] / W[i];
        }

        return 0;
      }

      template <typename ScalarT, typename IdxT>
      template <typename V>
      __attribute__((always_inline)) inline void
      SkinEffect<ScalarT, IdxT>::evaluateResidualKernel(V* y,
                                                        V*,
                                                        V*,
                                                        V* f) const
      {
        using namespace GridKit::LinearAlgebra;

        const Layout& L     = layout_;
        const auto    omega = static_cast<ScalarT>(coordinate());
        auto          Q     = slice(y, L.q, L.K, L.S);
        auto          G     = slice(y, L.g, L.K);
        auto          H     = slice(y, L.h, L.K);
        auto          W     = slice(y, L.w, L.K);
        auto          R     = slice(y, L.r, L.K);
        auto          Li    = slice(y, L.l, L.K);
        auto          QF    = slice(f, L.q, L.K, L.S);
        auto          GF    = slice(f, L.g, L.K);
        auto          HF    = slice(f, L.h, L.K);
        auto          WF    = slice(f, L.w, L.K);
        auto          RF    = slice(f, L.r, L.K);
        auto          LF    = slice(f, L.l, L.K);
        auto          Qref  = matrix(L.K,
                           L.S,
                           [&](IdxT i, IdxT k)
                           {
                             const auto rb = branchResistance(i, k);
                             const auto lb = branchInductance(i);
                             return rb * rb + omega * omega * lb * lb;
                           });
        auto          Gref  = vector(L.K,
                           [&](IdxT i)
                           {
                             V value = 0.0;
                             for (IdxT k = 0; k < L.S; ++k)
                             {
                               value += branchResistance(i, k) / eval(Q, i, k);
                             }
                             return value;
                           });
        auto          Href  = vector(L.K,
                           [&](IdxT i)
                           {
                             V          value = 0.0;
                             const auto lb    = branchInductance(i);
                             for (IdxT k = 0; k < L.S; ++k)
                             {
                               value += lb / eval(Q, i, k);
                             }
                             return value;
                           });

        equation(QF) = -Q + Qref;
        equation(GF) = -G + Gref;
        equation(HF) = -H + Href;
        equation(WF) = -W + hadamard(G, G) + omega * omega * hadamard(H, H);
        equation(RF) = -hadamard(W, R) + G;
        equation(LF) = -hadamard(W, Li) + H;
      }

      template <typename ScalarT, typename IdxT>
      void SkinEffect<ScalarT, IdxT>::evaluateInternalResidual(ScalarT* y,
                                                               ScalarT* yp,
                                                               ScalarT* u,
                                                               ScalarT* f) const
      {
        evaluateResidualKernel<ScalarT>(y, yp, u, f);
      }

      template <typename ScalarT, typename IdxT>
      int SkinEffect<ScalarT, IdxT>::evaluateResidual()
      {
        evaluateInternalResidual(y_, yp_, nullptr, f_);
        return 0;
      }

      template <typename ScalarT, typename IdxT>
      int SkinEffect<ScalarT, IdxT>::evaluateJacobian()
      {
        return this->template evaluateLocalAlgebraicJacobian<SkinEffect<ScalarT, IdxT>>();
      }

      template class SkinEffect<double, long int>;
      template class SkinEffect<double, size_t>;
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
