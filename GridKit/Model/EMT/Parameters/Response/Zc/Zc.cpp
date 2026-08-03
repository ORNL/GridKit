/**
 * @file Zc.cpp
 *
 * @brief Initialization, equation residual, and sparse-jet Jacobian for the Zc model.
 *
 */

#include "Zc.hpp"

#include <complex>
#include <stdexcept>

#include <GridKit/LinearAlgebra/DenseExpression.hpp>

#include <Eigen/LU>
#include <unsupported/Eigen/MatrixFunctions>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename ScalarT, typename IdxT>
      Zc<ScalarT, IdxT>::Zc(const SeriesImpedance<ScalarT, IdxT>& series,
                            const ShuntAdmittance<ScalarT, IdxT>& shunt)
        : series_(series),
          shunt_(shunt),
          layout_(series.R().rows)
      {
      }

      template <typename ScalarT, typename IdxT>
      int Zc<ScalarT, IdxT>::initialize()
      {
        using ComplexT = std::complex<ScalarT>;
        using MatrixXc = Eigen::Matrix<ComplexT, Eigen::Dynamic, Eigen::Dynamic>;

        auto eigenIndex = [](IdxT i)
        { return static_cast<Eigen::Index>(i); };

        const Layout& layout = layout_;
        const ScalarT omega  = static_cast<ScalarT>(coordinate());
        auto          Rc     = state(layout.rc, layout.K, layout.K);
        auto          Xc     = state(layout.xc, layout.K, layout.K);
        auto          R      = input(series_.R());
        auto          L      = input(series_.L());
        auto          G      = input(shunt_.G());
        auto          C      = input(shunt_.C());

        MatrixXc Z(eigenIndex(layout.K), eigenIndex(layout.K));
        MatrixXc Y(eigenIndex(layout.K), eigenIndex(layout.K));
        for (IdxT i = 0; i < layout.K; ++i)
        {
          for (IdxT j = 0; j < layout.K; ++j)
          {
            Z(eigenIndex(i), eigenIndex(j)) = ComplexT{R(i, j), omega * L(i, j)};
            Y(eigenIndex(i), eigenIndex(j)) = ComplexT{G(i, j), omega * C(i, j)};
          }
        }

        Eigen::FullPivLU<MatrixXc> y_lu(Y);
        if (!y_lu.isInvertible())
        {
          throw std::runtime_error("Zc initialization requires an invertible shunt admittance matrix");
        }

        const MatrixXc quotient = y_lu.solve(Z);
        const MatrixXc root     = quotient.sqrt();
        for (IdxT i = 0; i < layout.K; ++i)
        {
          for (IdxT j = 0; j < layout.K; ++j)
          {
            const auto value = root(eigenIndex(i), eigenIndex(j));
            Rc(i, j)         = value.real();
            Xc(i, j)         = value.imag();
          }
        }

        return 0;
      }

      template <typename ScalarT, typename IdxT>
      template <typename V>
      __attribute__((always_inline)) inline void Zc<ScalarT, IdxT>::evaluateResidualKernel(V* y,
                                                                                           V*,
                                                                                           V* u,
                                                                                           V* f) const
      {
        using namespace GridKit::LinearAlgebra;

        const Layout& layout = layout_;
        const auto    omega  = static_cast<ScalarT>(coordinate());
        auto          Rc     = slice(y, layout.rc, layout.K, layout.K);
        auto          Xc     = slice(y, layout.xc, layout.K, layout.K);
        auto          RFc    = slice(f, layout.rc, layout.K, layout.K);
        auto          XFc    = slice(f, layout.xc, layout.K, layout.K);
        auto          R      = input(u, series_.R());
        auto          L      = input(u, series_.L());
        auto          G      = input(u, shunt_.G());
        auto          C      = input(u, shunt_.C());

        auto X   = omega * L;
        auto BC  = omega * C;
        auto Z2r = Rc * Rc - Xc * Xc;
        auto Z2i = Rc * Xc + Xc * Rc;

        equation(RFc) = -R + G * Z2r - BC * Z2i;
        equation(XFc) = -X + G * Z2i + BC * Z2r;
      }

      template <typename ScalarT, typename IdxT>
      void Zc<ScalarT, IdxT>::evaluateInternalResidual(ScalarT* y,
                                                       ScalarT* yp,
                                                       ScalarT* u,
                                                       ScalarT* f) const
      {
        evaluateResidualKernel<ScalarT>(y, yp, u, f);
      }

      template <typename ScalarT, typename IdxT>
      int Zc<ScalarT, IdxT>::evaluateResidual()
      {
        evaluateInternalResidual(y_, yp_, nullptr, f_);
        return 0;
      }

      template <typename ScalarT, typename IdxT>
      int Zc<ScalarT, IdxT>::evaluateJacobian()
      {
        return this->template evaluateInputAlgebraicJacobian<Zc<ScalarT, IdxT>>();
      }

      template class Zc<double, long int>;
      template class Zc<double, size_t>;
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
