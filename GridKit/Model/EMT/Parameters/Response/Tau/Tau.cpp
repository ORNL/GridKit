/**
 * @file Tau.cpp
 *
 * @brief Initialization, equation residual, and sparse-jet Jacobian for the Tau model.
 *
 */

#include "Tau.hpp"

#include <GridKit/LinearAlgebra/DenseExpression.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename ScalarT, typename IdxT>
      Tau<ScalarT, IdxT>::Tau(const Gamma<ScalarT, IdxT>& gamma,
                              const Path<ScalarT, IdxT>&  path)
        : gamma_(gamma),
          path_(path),
          layout_(gamma.betaM().rows)
      {
      }

      template <typename ScalarT, typename IdxT>
      int Tau<ScalarT, IdxT>::initialize()
      {
        const ScalarT omega = static_cast<ScalarT>(coordinate());
        const ScalarT ell   = path_.length();
        auto          tau   = state(layout_.tau, layout_.K);
        auto          beta  = input(gamma_.betaM());

        for (IdxT m = 0; m < layout_.K; ++m)
        {
          tau[m] = ell * beta[m] / omega;
        }

        return 0;
      }

      template <typename ScalarT, typename IdxT>
      template <typename V>
      __attribute__((always_inline)) inline void
      Tau<ScalarT, IdxT>::evaluateResidualKernel(V* y,
                                                 V*,
                                                 V* u,
                                                 V* f) const
      {
        using namespace GridKit::LinearAlgebra;

        const auto omega = static_cast<ScalarT>(coordinate());
        const auto ell   = path_.length();
        auto       tau   = slice(y, layout_.tau, layout_.K);
        auto       Fres  = slice(f, layout_.tau, layout_.K);
        auto       beta  = input(u, gamma_.betaM());

        equation(Fres) = -tau + (ell / omega) * beta;
      }

      template <typename ScalarT, typename IdxT>
      void Tau<ScalarT, IdxT>::evaluateInternalResidual(ScalarT* y,
                                                        ScalarT* yp,
                                                        ScalarT* u,
                                                        ScalarT* f) const
      {
        evaluateResidualKernel<ScalarT>(y, yp, u, f);
      }

      template <typename ScalarT, typename IdxT>
      int Tau<ScalarT, IdxT>::evaluateResidual()
      {
        evaluateInternalResidual(y_, yp_, nullptr, f_);
        return 0;
      }

      template <typename ScalarT, typename IdxT>
      int Tau<ScalarT, IdxT>::evaluateJacobian()
      {
        return this->template evaluateInputAlgebraicJacobian<Tau<ScalarT, IdxT>>();
      }

      template class Tau<double, long int>;
      template class Tau<double, size_t>;
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
