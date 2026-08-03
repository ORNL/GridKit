/**
 * @file Yc.cpp
 *
 * @brief Initialization, equation residual, and sparse-jet Jacobian for the Yc model.
 *
 */

#include "Yc.hpp"

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
      Yc<ScalarT, IdxT>::Yc(const SeriesView<ScalarT, IdxT>& series,
                            const ShuntView<ScalarT, IdxT>&  shunt)
        : series_(series),
          shunt_(shunt),
          layout_(series.R().rows)
      {
      }

      template <typename ScalarT, typename IdxT>
      int Yc<ScalarT, IdxT>::initialize()
      {
        using ComplexT = std::complex<ScalarT>;
        using MatrixXc = Eigen::Matrix<ComplexT, Eigen::Dynamic, Eigen::Dynamic>;

        auto eigenIndex = [](IdxT i)
        { return static_cast<Eigen::Index>(i); };

        const Layout& layout = layout_;
        const ScalarT omega  = static_cast<ScalarT>(coordinate());
        auto          Gc     = state(layout.gc, layout.K, layout.K);
        auto          Bc     = state(layout.bc, layout.K, layout.K);
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

        Eigen::FullPivLU<MatrixXc> z_lu(Z);
        if (!z_lu.isInvertible())
        {
          throw std::runtime_error("Yc initialization requires an invertible series impedance matrix");
        }

        // Yc = Z^-1 sqrt(Z Y) solves the sandwich equation Yc Z Yc = Y
        // exactly, so the residual vanishes at this point and the result
        // is symmetric whenever Z and Y are.
        const MatrixXc product    = Z * Y;
        const MatrixXc admittance = z_lu.solve(product.sqrt());
        for (IdxT i = 0; i < layout.K; ++i)
        {
          for (IdxT j = 0; j < layout.K; ++j)
          {
            const auto value = admittance(eigenIndex(i), eigenIndex(j));
            Gc(i, j)         = value.real();
            Bc(i, j)         = value.imag();
          }
        }

        return 0;
      }

      template <typename ScalarT, typename IdxT>
      template <typename V>
      __attribute__((always_inline)) inline void Yc<ScalarT, IdxT>::evaluateResidualKernel(V* y,
                                                                                           V*,
                                                                                           V* u,
                                                                                           V* f) const
      {
        using namespace GridKit::LinearAlgebra;

        const Layout& layout = layout_;
        const auto    omega  = static_cast<ScalarT>(coordinate());
        auto          Gc     = slice(y, layout.gc, layout.K, layout.K);
        auto          Bc     = slice(y, layout.bc, layout.K, layout.K);
        auto          GFc    = slice(f, layout.gc, layout.K, layout.K);
        auto          BFc    = slice(f, layout.bc, layout.K, layout.K);
        auto          R      = input(u, series_.R());
        auto          L      = input(u, series_.L());
        auto          G      = input(u, shunt_.G());
        auto          C      = input(u, shunt_.C());

        auto B  = omega * C;
        auto XL = omega * L;

        // The sandwich equation 0 = -Y + Yc Z Yc, split over P = Z Yc.
        // The commuting form Z Yc Yc holds only when Z and Y commute and
        // loses the exact symmetry of the characteristic admittance.
        auto Pr = R * Gc - XL * Bc;
        auto Pi = R * Bc + XL * Gc;

        equation(GFc) = -G + Gc * Pr - Bc * Pi;
        equation(BFc) = -B + Gc * Pi + Bc * Pr;
      }

      template <typename ScalarT, typename IdxT>
      void Yc<ScalarT, IdxT>::evaluateInternalResidual(ScalarT* y,
                                                       ScalarT* yp,
                                                       ScalarT* u,
                                                       ScalarT* f) const
      {
        evaluateResidualKernel<ScalarT>(y, yp, u, f);
      }

      template <typename ScalarT, typename IdxT>
      int Yc<ScalarT, IdxT>::evaluateResidual()
      {
        evaluateInternalResidual(y_, yp_, nullptr, f_);
        return 0;
      }

      template <typename ScalarT, typename IdxT>
      int Yc<ScalarT, IdxT>::evaluateJacobian()
      {
        return this->template evaluateInputAlgebraicJacobian<Yc<ScalarT, IdxT>>();
      }

      template class Yc<double, long int>;
      template class Yc<double, size_t>;
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
