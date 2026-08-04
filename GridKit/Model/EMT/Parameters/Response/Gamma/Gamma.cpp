/**
 * @file Gamma.cpp
 *
 * @brief Initialization, equation residual, and sparse-jet Jacobian for the Gamma model.
 *
 */

#include "Gamma.hpp"

#include <complex>

#include <GridKit/LinearAlgebra/DenseExpression.hpp>

#include <unsupported/Eigen/MatrixFunctions>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename ScalarT, typename IdxT>
      Gamma<ScalarT, IdxT>::Gamma(const SeriesView<ScalarT, IdxT>& series,
                                  const ShuntView<ScalarT, IdxT>&  shunt)
        : series_(series),
          shunt_(shunt),
          layout_(series.R().rows)
      {
      }

      /**
       * @brief Seed the principal square root of ZY; the residual then
       *        holds the states on this branch across the sweep.
       */
      template <typename ScalarT, typename IdxT>
      int Gamma<ScalarT, IdxT>::initialize()
      {
        using ComplexT = std::complex<ScalarT>;
        using MatrixXc = Eigen::Matrix<ComplexT, Eigen::Dynamic, Eigen::Dynamic>;

        auto eigenIndex = [](IdxT i)
        { return static_cast<Eigen::Index>(i); };

        const Layout& layout = layout_;
        const ScalarT omega  = static_cast<ScalarT>(coordinate());
        auto          Alpha  = state(layout.alpha, layout.K, layout.K);
        auto          Beta   = state(layout.beta, layout.K, layout.K);
        auto          R      = input(series_.R());
        auto          L      = input(series_.L());
        auto          G      = input(shunt_.G());
        auto          C      = input(shunt_.C());

        MatrixXc product(eigenIndex(layout.K), eigenIndex(layout.K));
        for (IdxT i = 0; i < layout.K; ++i)
        {
          for (IdxT j = 0; j < layout.K; ++j)
          {
            product(eigenIndex(i), eigenIndex(j)) = ComplexT{0.0, 0.0};
            for (IdxT k = 0; k < layout.K; ++k)
            {
              product(eigenIndex(i), eigenIndex(j)) +=
                  ComplexT{R(i, k), omega * L(i, k)}
                  * ComplexT{G(k, j), omega * C(k, j)};
            }
          }
        }

        const MatrixXc gamma = product.sqrt();
        for (IdxT i = 0; i < layout.K; ++i)
        {
          for (IdxT j = 0; j < layout.K; ++j)
          {
            Alpha(i, j) = gamma(eigenIndex(i), eigenIndex(j)).real();
            Beta(i, j)  = gamma(eigenIndex(i), eigenIndex(j)).imag();
          }
        }

        return 0;
      }

      template <typename ScalarT, typename IdxT>
      template <typename V>
      __attribute__((always_inline)) inline void
      Gamma<ScalarT, IdxT>::evaluateResidualKernel(V* y,
                                                   V*,
                                                   V* u,
                                                   V* f) const
      {
        using namespace GridKit::LinearAlgebra;

        const Layout& layout = layout_;
        const auto    omega  = static_cast<ScalarT>(coordinate());
        auto          Alpha  = slice(y, layout.alpha, layout.K, layout.K);
        auto          Beta   = slice(y, layout.beta, layout.K, layout.K);
        auto          Ares   = slice(f, layout.alpha, layout.K, layout.K);
        auto          Bres   = slice(f, layout.beta, layout.K, layout.K);
        auto          R      = input(u, series_.R());
        auto          L      = input(u, series_.L());
        auto          G      = input(u, shunt_.G());
        auto          C      = input(u, shunt_.C());

        auto XL = omega * L;
        auto BC = omega * C;

        equation(Ares) = Alpha * Alpha - Beta * Beta - R * G + XL * BC;
        equation(Bres) = Alpha * Beta + Beta * Alpha - R * BC - XL * G;
      }

      template <typename ScalarT, typename IdxT>
      void Gamma<ScalarT, IdxT>::evaluateInternalResidual(ScalarT* y,
                                                          ScalarT* yp,
                                                          ScalarT* u,
                                                          ScalarT* f) const
      {
        evaluateResidualKernel<ScalarT>(y, yp, u, f);
      }

      template <typename ScalarT, typename IdxT>
      void Gamma<ScalarT, IdxT>::tagDifferentiable(std::vector<bool>& tag) const
      {
        for (IdxT i = 0; i < size(); ++i)
        {
          tag[static_cast<size_t>(base_ + i)] = true;
        }
      }

      template <typename ScalarT, typename IdxT>
      int Gamma<ScalarT, IdxT>::evaluateResidual()
      {
        evaluateInternalResidual(y_, yp_, nullptr, f_);
        return 0;
      }

      template <typename ScalarT, typename IdxT>
      int Gamma<ScalarT, IdxT>::evaluateJacobian()
      {
        return this->template evaluateInputAlgebraicJacobian<Gamma<ScalarT, IdxT>>();
      }

      template class Gamma<double, long int>;
      template class Gamma<double, size_t>;
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
