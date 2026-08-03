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

        // Zc = Y^-1 sqrt(Y Z) solves the sandwich equation Zc Y Zc = Z
        // exactly, so the residual vanishes at this point and the result
        // is symmetric whenever Z and Y are.
        const MatrixXc product   = Y * Z;
        const MatrixXc impedance = y_lu.solve(product.sqrt());
        for (IdxT i = 0; i < layout.K; ++i)
        {
          for (IdxT j = 0; j < layout.K; ++j)
          {
            const auto value = impedance(eigenIndex(i), eigenIndex(j));
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

        auto X  = omega * L;
        auto BC = omega * C;

        // The sandwich equation 0 = -Z + Zc Y Zc, split over P = Y Zc.
        // The commuting form Y Zc Zc holds only when Z and Y commute and
        // loses the exact symmetry of the characteristic impedance.
        auto Pr = G * Rc - BC * Xc;
        auto Pi = G * Xc + BC * Rc;

        equation(RFc) = -R + Rc * Pr - Xc * Pi;
        equation(XFc) = -X + Rc * Pi + Xc * Pr;
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
