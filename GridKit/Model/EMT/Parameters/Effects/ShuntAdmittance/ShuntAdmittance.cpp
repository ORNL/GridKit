/**
 * @file ShuntAdmittance.cpp
 *
 * @brief Initialization, equation residual, and sparse-jet Jacobian for the ShuntAdmittance
 * model.
 *
 */

#include "ShuntAdmittance.hpp"

#include <GridKit/LinearAlgebra/DenseExpression.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename ScalarT, typename IdxT>
      ShuntAdmittance<ScalarT, IdxT>::ShuntAdmittance(
          const ShuntPotential<ScalarT, IdxT>& shunt)
        : shunt_(shunt),
          layout_(shunt.Gpot().rows)
      {
      }

      template <typename ScalarT, typename IdxT>
      int ShuntAdmittance<ScalarT, IdxT>::initialize()
      {
        const Layout& layout = layout_;
        auto          Gext   = state(layout.g_ext, layout.K, layout.K);
        auto          Cext   = state(layout.c_ext, layout.K, layout.K);
        auto          Gshunt = state(layout.g, layout.K, layout.K);
        auto          Cshunt = state(layout.c, layout.K, layout.K);
        auto          Gpot   = input(shunt_.Gpot());
        auto          Cpot   = input(shunt_.Cpot());

        for (IdxT i = 0; i < layout.K; ++i)
        {
          for (IdxT j = 0; j < layout.K; ++j)
          {
            Gext(i, j)   = 0.0;
            Cext(i, j)   = 0.0;
            Gshunt(i, j) = Gpot(i, j);
            Cshunt(i, j) = Cpot(i, j);
          }
        }

        return 0;
      }

      template <typename ScalarT, typename IdxT>
      template <typename V>
      __attribute__((always_inline)) inline void
      ShuntAdmittance<ScalarT, IdxT>::evaluateResidualKernel(V* y,
                                                             V*,
                                                             V* u,
                                                             V* f) const
      {
        using namespace GridKit::LinearAlgebra;

        const Layout& layout = layout_;
        auto          Gext   = slice(y, layout.g_ext, layout.K, layout.K);
        auto          Cext   = slice(y, layout.c_ext, layout.K, layout.K);
        auto          Gshunt = slice(y, layout.g, layout.K, layout.K);
        auto          Cshunt = slice(y, layout.c, layout.K, layout.K);
        auto          GFext  = slice(f, layout.g_ext, layout.K, layout.K);
        auto          CFext  = slice(f, layout.c_ext, layout.K, layout.K);
        auto          GF     = slice(f, layout.g, layout.K, layout.K);
        auto          CF     = slice(f, layout.c, layout.K, layout.K);
        auto          Gpot   = input(u, shunt_.Gpot());
        auto          Cpot   = input(u, shunt_.Cpot());

        equation(GFext) = -Gext;
        equation(CFext) = -Cext;
        equation(GF)    = -Gshunt + Gpot + Gext;
        equation(CF)    = -Cshunt + Cpot + Cext;
      }

      template <typename ScalarT, typename IdxT>
      void ShuntAdmittance<ScalarT, IdxT>::evaluateInternalResidual(ScalarT* y,
                                                                    ScalarT* yp,
                                                                    ScalarT* u,
                                                                    ScalarT* f) const
      {
        evaluateResidualKernel<ScalarT>(y, yp, u, f);
      }

      template <typename ScalarT, typename IdxT>
      int ShuntAdmittance<ScalarT, IdxT>::evaluateResidual()
      {
        evaluateInternalResidual(y_, yp_, nullptr, f_);
        return 0;
      }

      template <typename ScalarT, typename IdxT>
      int ShuntAdmittance<ScalarT, IdxT>::evaluateJacobian()
      {
        return this->template evaluateInputAlgebraicJacobian<ShuntAdmittance<ScalarT, IdxT>>();
      }

      template class ShuntAdmittance<double, long int>;
      template class ShuntAdmittance<double, size_t>;
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
