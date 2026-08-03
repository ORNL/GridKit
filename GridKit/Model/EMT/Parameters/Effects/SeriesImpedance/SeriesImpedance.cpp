/**
 * @file SeriesImpedance.cpp
 *
 * @brief Initialization, equation residual, and sparse-jet Jacobian for the SeriesImpedance
 * model.
 *
 */

#include "SeriesImpedance.hpp"

#include <GridKit/LinearAlgebra/DenseExpression.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename ScalarT, typename IdxT>
      SeriesImpedance<ScalarT, IdxT>::SeriesImpedance(
          const SkinEffect<ScalarT, IdxT>&          skin,
          const GeometricInductance<ScalarT, IdxT>& geometric,
          const Carson<ScalarT, IdxT>&              carson)
        : skin_(skin),
          geometric_(geometric),
          carson_(carson),
          layout_(skin.rskin().rows)
      {
      }

      template <typename ScalarT, typename IdxT>
      int SeriesImpedance<ScalarT, IdxT>::initialize()
      {
        const Layout& layout  = layout_;
        auto          Rint    = state(layout.r_int, layout.K, layout.K);
        auto          Lint    = state(layout.l_int, layout.K, layout.K);
        auto          Rext    = state(layout.r_ext, layout.K, layout.K);
        auto          Lext    = state(layout.l_ext, layout.K, layout.K);
        auto          Rseries = state(layout.r, layout.K, layout.K);
        auto          Lseries = state(layout.l, layout.K, layout.K);
        auto          rskin   = input(skin_.rskin());
        auto          lskin   = input(skin_.lskin());
        auto          Lgeo    = input(geometric_.Lgeo());
        auto          Rcarson = input(carson_.Rcarson());
        auto          Lcarson = input(carson_.Lcarson());

        for (IdxT i = 0; i < layout.K; ++i)
        {
          for (IdxT j = 0; j < layout.K; ++j)
          {
            Rint(i, j)    = (i == j) ? rskin[i] : 0.0;
            Lint(i, j)    = (i == j) ? lskin[i] : 0.0;
            Rext(i, j)    = Rcarson(i, j);
            Lext(i, j)    = Lgeo(i, j) + Lcarson(i, j);
            Rseries(i, j) = Rint(i, j) + Rext(i, j);
            Lseries(i, j) = Lint(i, j) + Lext(i, j);
          }
        }

        return 0;
      }

      template <typename ScalarT, typename IdxT>
      template <typename V>
      __attribute__((always_inline)) inline void
      SeriesImpedance<ScalarT, IdxT>::evaluateResidualKernel(V* y,
                                                             V*,
                                                             V* u,
                                                             V* f) const
      {
        using namespace GridKit::LinearAlgebra;

        const Layout& layout  = layout_;
        auto          Rint    = slice(y, layout.r_int, layout.K, layout.K);
        auto          Lint    = slice(y, layout.l_int, layout.K, layout.K);
        auto          Rext    = slice(y, layout.r_ext, layout.K, layout.K);
        auto          Lext    = slice(y, layout.l_ext, layout.K, layout.K);
        auto          Rseries = slice(y, layout.r, layout.K, layout.K);
        auto          Lseries = slice(y, layout.l, layout.K, layout.K);
        auto          RFint   = slice(f, layout.r_int, layout.K, layout.K);
        auto          LFint   = slice(f, layout.l_int, layout.K, layout.K);
        auto          RFext   = slice(f, layout.r_ext, layout.K, layout.K);
        auto          LFext   = slice(f, layout.l_ext, layout.K, layout.K);
        auto          RF      = slice(f, layout.r, layout.K, layout.K);
        auto          LF      = slice(f, layout.l, layout.K, layout.K);
        auto          rskin   = input(u, skin_.rskin());
        auto          lskin   = input(u, skin_.lskin());
        auto          Lgeo    = input(u, geometric_.Lgeo());
        auto          Rcarson = input(u, carson_.Rcarson());
        auto          Lcarson = input(u, carson_.Lcarson());

        equation(RFint) = -Rint + diag(rskin);
        equation(LFint) = -Lint + diag(lskin);
        equation(RFext) = -Rext + Rcarson;
        equation(LFext) = -Lext + Lgeo + Lcarson;
        equation(RF)    = -Rseries + Rint + Rext;
        equation(LF)    = -Lseries + Lint + Lext;
      }

      template <typename ScalarT, typename IdxT>
      void SeriesImpedance<ScalarT, IdxT>::evaluateInternalResidual(ScalarT* y,
                                                                    ScalarT* yp,
                                                                    ScalarT* u,
                                                                    ScalarT* f) const
      {
        evaluateResidualKernel<ScalarT>(y, yp, u, f);
      }

      template <typename ScalarT, typename IdxT>
      int SeriesImpedance<ScalarT, IdxT>::evaluateResidual()
      {
        evaluateInternalResidual(y_, yp_, nullptr, f_);
        return 0;
      }

      template <typename ScalarT, typename IdxT>
      int SeriesImpedance<ScalarT, IdxT>::evaluateJacobian()
      {
        return this->template evaluateInputAlgebraicJacobian<SeriesImpedance<ScalarT, IdxT>>();
      }

      template class SeriesImpedance<double, long int>;
      template class SeriesImpedance<double, size_t>;
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
