/**
 * @file Carson.cpp
 *
 * @brief Initialization, equation residual, and sparse-jet Jacobian for the Carson model.
 *
 */

#include "Carson.hpp"

#include <cmath>
#include <stdexcept>

#include <GridKit/LinearAlgebra/DenseExpression.hpp>
#include <GridKit/Model/EMT/Constants.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename ScalarT, typename IdxT>
      Carson<ScalarT, IdxT>::Carson(const Tower<ScalarT, IdxT>& tower,
                                    CarsonData<ScalarT, IdxT>   data)
        : tower_(tower),
          data_(std::move(data)),
          layout_(tower.conductorCount())
      {
      }

      template <typename ScalarT, typename IdxT>
      template <typename V>
      typename Carson<ScalarT, IdxT>::template PairT<V>
      Carson<ScalarT, IdxT>::evaluatePair(IdxT i, IdxT j, V gr, V gi) const
      {
        using std::atan2;
        using std::log;

        const auto omega = static_cast<ScalarT>(coordinate());
        const auto mu0   = Constants::mu0<ScalarT>();
        const auto pi    = Constants::pi<ScalarT>();
        const auto d     = tower_.horizontalDistance(i, j);
        const auto s     = tower_.height(i) + tower_.height(j);
        const auto dp    = tower_.imageDistance(i, j);
        const auto gg    = gr * gr + gi * gi;

        const auto p_real = gr / gg;
        const auto p_imag = -gi / gg;
        const auto h_real = s + 2.0 * p_real;
        const auto h_imag = 2.0 * p_imag;
        const auto A      = d * d + h_real * h_real - h_imag * h_imag;
        const auto B      = 2.0 * h_real * h_imag;
        const auto M      = A * A + B * B;

        const auto theta  = atan2(B, A);
        const auto lambda = 0.25 * log(M) - log(dp);

        PairT<V> pair;
        pair.r = -omega * mu0 * theta / (4.0 * pi);
        pair.l = mu0 * lambda / (2.0 * pi);
        return pair;
      }

      template <typename ScalarT, typename IdxT>
      int Carson<ScalarT, IdxT>::initialize()
      {
        const Layout& L     = layout_;
        const ScalarT omega = static_cast<ScalarT>(coordinate());
        if (!std::isfinite(static_cast<double>(sigma())) || sigma() < static_cast<ScalarT>(0.0))
        {
          throw std::runtime_error("Carson earth conductivity must be nonnegative and finite");
        }
        if (!std::isfinite(static_cast<double>(epsilon()))
            || epsilon() < Constants::epsilon0<ScalarT>())
        {
          throw std::runtime_error("Carson earth permittivity must be at least vacuum permittivity");
        }
        if (!std::isfinite(static_cast<double>(omega)) || omega <= static_cast<ScalarT>(0.0))
        {
          throw std::runtime_error("Carson angular frequency must be positive and finite");
        }

        const ScalarT a  = omega * omega * Constants::mu0<ScalarT>() * epsilon();
        const ScalarT b  = omega * Constants::mu0<ScalarT>() * sigma();
        auto          R  = state(L.r, L.K, L.K);
        auto          Li = state(L.l, L.K, L.K);

        y_[L.eta]     = std::sqrt(a * a + b * b);
        y_[L.gamma_r] = std::sqrt(0.5 * (y_[L.eta] - a));
        y_[L.gamma_i] = std::sqrt(0.5 * (y_[L.eta] + a));

        for (IdxT i = 0; i < L.K; ++i)
        {
          for (IdxT j = 0; j < L.K; ++j)
          {
            const Pair pair = evaluatePair(i, j, y_[L.gamma_r], y_[L.gamma_i]);
            R(i, j)         = pair.r;
            Li(i, j)        = pair.l;
          }
        }

        return 0;
      }

      template <typename ScalarT, typename IdxT>
      template <typename V>
      __attribute__((always_inline)) inline void Carson<ScalarT, IdxT>::evaluateResidualKernel(V* y,
                                                                                               V*,
                                                                                               V*,
                                                                                               V* f) const
      {
        using namespace GridKit::LinearAlgebra;

        const Layout& L       = layout_;
        const auto    omega   = static_cast<ScalarT>(coordinate());
        const auto    mu0     = Constants::mu0<ScalarT>();
        const auto    pi      = Constants::pi<ScalarT>();
        const auto    a       = omega * omega * mu0 * epsilon();
        const auto    b       = omega * mu0 * sigma();
        auto          Eta     = slice(y, L.eta, 1);
        auto          gammaR  = slice(y, L.gamma_r, 1);
        auto          gammaI  = slice(y, L.gamma_i, 1);
        auto          R       = slice(y, L.r, L.K, L.K);
        auto          Li      = slice(y, L.l, L.K, L.K);
        auto          etaF    = slice(f, L.eta, 1);
        auto          gammaRF = slice(f, L.gamma_r, 1);
        auto          gammaIF = slice(f, L.gamma_i, 1);
        auto          RF      = slice(f, L.r, L.K, L.K);
        auto          LF      = slice(f, L.l, L.K, L.K);

        const auto gr     = eval(gammaR, 0);
        const auto gi     = eval(gammaI, 0);
        auto       etaRef = vector(1,
                             [&](IdxT)
                             {
                               return a * a + b * b;
                             });
        auto       aRef   = vector(1,
                           [&](IdxT)
                           {
                             return a;
                           });
        auto       Rref   = matrix(L.K,
                           L.K,
                           [&](IdxT i, IdxT j)
                           {
                             return evaluatePair(i, j, gr, gi).r;
                           });
        auto       Lref   = matrix(L.K,
                           L.K,
                           [&](IdxT i, IdxT j)
                           {
                             return evaluatePair(i, j, gr, gi).l;
                           });

        equation(etaF)    = -hadamard(Eta, Eta) + etaRef;
        equation(gammaRF) = -2.0 * hadamard(gammaR, gammaR) + Eta - aRef;
        equation(gammaIF) = -2.0 * hadamard(gammaI, gammaI) + Eta + aRef;
        equation(RF)      = -4.0 * pi * R + 4.0 * pi * Rref;
        equation(LF)      = -2.0 * pi * Li + 2.0 * pi * Lref;
      }

      template <typename ScalarT, typename IdxT>
      void Carson<ScalarT, IdxT>::evaluateInternalResidual(ScalarT* y,
                                                           ScalarT* yp,
                                                           ScalarT* u,
                                                           ScalarT* f) const
      {
        evaluateResidualKernel<ScalarT>(y, yp, u, f);
      }

      template <typename ScalarT, typename IdxT>
      int Carson<ScalarT, IdxT>::evaluateResidual()
      {
        evaluateInternalResidual(y_, yp_, nullptr, f_);
        return 0;
      }

      template <typename ScalarT, typename IdxT>
      int Carson<ScalarT, IdxT>::evaluateJacobian()
      {
        return this->template evaluateLocalAlgebraicJacobian<Carson<ScalarT, IdxT>>();
      }

      template class Carson<double, long int>;
      template class Carson<double, size_t>;
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
