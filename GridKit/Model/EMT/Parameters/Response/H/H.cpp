/**
 * @file H.cpp
 *
 * @brief Initialization, equation residual, and sparse-jet Jacobian for the H model.
 *
 */

#include "H.hpp"

#include <cmath>

#include <GridKit/LinearAlgebra/DenseExpression.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename ScalarT, typename IdxT>
      H<ScalarT, IdxT>::H(const Gamma<ScalarT, IdxT>& gamma,
                          const Path<ScalarT, IdxT>&  path)
        : gamma_(gamma),
          path_(path),
          layout_(gamma.alphaM().rows)
      {
      }

      template <typename ScalarT, typename IdxT>
      int H<ScalarT, IdxT>::initialize()
      {
        auto Hr    = state(layout_.hr, layout_.K);
        auto Hi    = state(layout_.hi, layout_.K);
        auto Hrd   = derivative(layout_.hr, layout_.K);
        auto Hid   = derivative(layout_.hi, layout_.K);
        auto ell   = path_.length();
        auto alpha = input(gamma_.alphaM());
        auto beta  = input(gamma_.betaM());
        auto ad    = input(gamma_.alphaMDot());
        auto bd    = input(gamma_.betaMDot());

        for (IdxT m = 0; m < layout_.K; ++m)
        {
          const ScalarT attenuation = std::exp(-alpha[m] * ell);
          const ScalarT phase       = beta[m] * ell;
          Hr[m]                     = attenuation * std::cos(phase);
          Hi[m]                     = -attenuation * std::sin(phase);
          Hrd[m]                    = -ell * (ad[m] * Hr[m] - bd[m] * Hi[m]);
          Hid[m]                    = -ell * (ad[m] * Hi[m] + bd[m] * Hr[m]);
        }

        return 0;
      }

      template <typename ScalarT, typename IdxT>
      template <typename V>
      __attribute__((always_inline)) inline void H<ScalarT, IdxT>::evaluateResidualKernel(V* y,
                                                                                          V* yp,
                                                                                          V* u,
                                                                                          V* f) const
      {
        using namespace GridKit::LinearAlgebra;

        auto Hr  = slice(y, layout_.hr, layout_.K);
        auto Hi  = slice(y, layout_.hi, layout_.K);
        auto Hrd = slice(yp, layout_.hr, layout_.K);
        auto Hid = slice(yp, layout_.hi, layout_.K);
        auto Fhr = slice(f, layout_.hr, layout_.K);
        auto Fhi = slice(f, layout_.hi, layout_.K);
        auto ad  = input(u, gamma_.alphaMDot());
        auto bd  = input(u, gamma_.betaMDot());
        auto ell = path_.length();

        equation(Fhr) = -Hrd - ell * (hadamard(ad, Hr) - hadamard(bd, Hi));
        equation(Fhi) = -Hid - ell * (hadamard(ad, Hi) + hadamard(bd, Hr));
      }

      template <typename ScalarT, typename IdxT>
      void H<ScalarT, IdxT>::evaluateInternalResidual(ScalarT* y,
                                                      ScalarT* yp,
                                                      ScalarT* u,
                                                      ScalarT* f) const
      {
        evaluateResidualKernel<ScalarT>(y, yp, u, f);
      }

      template <typename ScalarT, typename IdxT>
      int H<ScalarT, IdxT>::evaluateResidual()
      {
        evaluateInternalResidual(y_, yp_, nullptr, f_);
        return 0;
      }

      template <typename ScalarT, typename IdxT>
      int H<ScalarT, IdxT>::evaluateJacobian()
      {
        return this->template evaluateInputDifferentialJacobian<H<ScalarT, IdxT>>();
      }

      template class H<double, long int>;
      template class H<double, size_t>;
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
