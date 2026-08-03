/**
 * @file Gamma.cpp
 *
 * @brief Initialization, equation residual, and sparse-jet Jacobian for the Gamma model.
 *
 */

#include "Gamma.hpp"

#include <algorithm>
#include <complex>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <GridKit/LinearAlgebra/DenseExpression.hpp>

#include <Eigen/Eigenvalues>

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

      template <typename ScalarT, typename IdxT>
      int Gamma<ScalarT, IdxT>::initialize()
      {
        using ComplexT = std::complex<ScalarT>;
        using MatrixXc = Eigen::Matrix<ComplexT, Eigen::Dynamic, Eigen::Dynamic>;
        using VectorXc = Eigen::Matrix<ComplexT, Eigen::Dynamic, 1>;

        auto eigenIndex = [](IdxT i)
        { return static_cast<Eigen::Index>(i); };

        const Layout& layout = layout_;
        const ScalarT omega  = static_cast<ScalarT>(coordinate());
        auto          Alpha  = state(layout.alpha, layout.K, layout.K);
        auto          Beta   = state(layout.beta, layout.K, layout.K);
        auto          am     = state(layout.a, layout.K);
        auto          bm     = state(layout.b, layout.K);
        auto          TvR    = state(layout.tv_r, layout.K, layout.K);
        auto          TvI    = state(layout.tv_i, layout.K, layout.K);
        auto          TiR    = state(layout.ti_r, layout.K, layout.K);
        auto          TiI    = state(layout.ti_i, layout.K, layout.K);
        auto          R      = input(series_.R());
        auto          L      = input(series_.L());
        auto          G      = input(shunt_.G());
        auto          C      = input(shunt_.C());

        MatrixXc gamma2(eigenIndex(layout.K), eigenIndex(layout.K));
        gamma2.setZero();
        for (IdxT i = 0; i < layout.K; ++i)
        {
          for (IdxT j = 0; j < layout.K; ++j)
          {
            for (IdxT k = 0; k < layout.K; ++k)
            {
              gamma2(eigenIndex(i), eigenIndex(j)) +=
                  ComplexT{R(i, k), omega * L(i, k)}
                  * ComplexT{G(k, j), omega * C(k, j)};
            }
          }
        }

        Eigen::ComplexEigenSolver<MatrixXc> solver(gamma2, true);
        if (solver.info() != Eigen::Success)
        {
          throw std::runtime_error("Gamma eigensolve failed");
        }

        auto     lambda2 = solver.eigenvalues();
        MatrixXc V       = solver.eigenvectors();

        std::vector<IdxT> order(static_cast<size_t>(layout.K));
        std::iota(order.begin(), order.end(), IdxT{0});
        std::sort(order.begin(),
                  order.end(),
                  [&](IdxT lhs, IdxT rhs)
                  {
                    const ComplexT l = std::sqrt(lambda2[eigenIndex(lhs)]);
                    const ComplexT r = std::sqrt(lambda2[eigenIndex(rhs)]);
                    return (l.imag() < r.imag()) || (l.imag() == r.imag() && l.real() < r.real());
                  });

        MatrixXc Vs(eigenIndex(layout.K), eigenIndex(layout.K));
        VectorXc lambda(eigenIndex(layout.K));
        for (IdxT m = 0; m < layout.K; ++m)
        {
          lambda[eigenIndex(m)] = std::sqrt(lambda2[eigenIndex(order[static_cast<size_t>(m)])]);
          Vs.col(eigenIndex(m)) = V.col(eigenIndex(order[static_cast<size_t>(m)]));
        }

        MatrixXc Wstar = Vs.inverse();

        for (IdxT m = 0; m < layout.K; ++m)
        {
          Eigen::Index pivot = 0;
          ScalarT      norm  = 0.0;
          for (IdxT i = 0; i < layout.K; ++i)
          {
            const ScalarT candidate = std::abs(Vs(eigenIndex(i), eigenIndex(m)));
            if (candidate > norm)
            {
              pivot = eigenIndex(i);
              norm  = candidate;
            }
          }

          if (norm > 0.0)
          {
            const ComplexT phase      = std::conj(Vs(pivot, eigenIndex(m))) / norm;
            Vs.col(eigenIndex(m))    *= phase;
            Wstar.row(eigenIndex(m)) /= phase;
          }
        }

        const MatrixXc gamma = Vs * lambda.asDiagonal() * Wstar;
        const MatrixXc Ti    = Wstar.adjoint();

        for (IdxT i = 0; i < layout.K; ++i)
        {
          am[i] = lambda[eigenIndex(i)].real();
          bm[i] = lambda[eigenIndex(i)].imag();

          for (IdxT j = 0; j < layout.K; ++j)
          {
            Alpha(i, j) = gamma(eigenIndex(i), eigenIndex(j)).real();
            Beta(i, j)  = gamma(eigenIndex(i), eigenIndex(j)).imag();
            TvR(i, j)   = Vs(eigenIndex(i), eigenIndex(j)).real();
            TvI(i, j)   = Vs(eigenIndex(i), eigenIndex(j)).imag();
            TiR(i, j)   = Ti(eigenIndex(i), eigenIndex(j)).real();
            TiI(i, j)   = Ti(eigenIndex(i), eigenIndex(j)).imag();
          }
        }

        return 0;
      }

      template <typename ScalarT, typename IdxT>
      template <typename V>
      __attribute__((always_inline)) inline void
      Gamma<ScalarT, IdxT>::evaluateResidualKernel(V* y,
                                                   V* yp,
                                                   V* u,
                                                   V* f) const
      {
        using namespace GridKit::LinearAlgebra;

        const Layout& layout = layout_;
        const auto    omega  = static_cast<ScalarT>(coordinate());
        auto          Alpha  = slice(y, layout.alpha, layout.K, layout.K);
        auto          Beta   = slice(y, layout.beta, layout.K, layout.K);
        auto          am     = slice(y, layout.a, layout.K);
        auto          bm     = slice(y, layout.b, layout.K);
        auto          TvR    = slice(y, layout.tv_r, layout.K, layout.K);
        auto          TvI    = slice(y, layout.tv_i, layout.K, layout.K);
        auto          TiR    = slice(y, layout.ti_r, layout.K, layout.K);
        auto          TiI    = slice(y, layout.ti_i, layout.K, layout.K);
        auto          Alphad = slice(yp, layout.alpha, layout.K, layout.K);
        auto          Betad  = slice(yp, layout.beta, layout.K, layout.K);
        auto          amd    = slice(yp, layout.a, layout.K);
        auto          bmd    = slice(yp, layout.b, layout.K);
        auto          TvRd   = slice(yp, layout.tv_r, layout.K, layout.K);
        auto          TvId   = slice(yp, layout.tv_i, layout.K, layout.K);
        auto          TiRd   = slice(yp, layout.ti_r, layout.K, layout.K);
        auto          TiId   = slice(yp, layout.ti_i, layout.K, layout.K);
        auto          Ares   = slice(f, layout.alpha, layout.K, layout.K);
        auto          Bres   = slice(f, layout.beta, layout.K, layout.K);
        auto          ares   = slice(f, layout.a, layout.K);
        auto          bres   = slice(f, layout.b, layout.K);
        auto          TvRres = slice(f, layout.tv_r, layout.K, layout.K);
        auto          TvIres = slice(f, layout.tv_i, layout.K, layout.K);
        auto          TiRres = slice(f, layout.ti_r, layout.K, layout.K);
        auto          TiIres = slice(f, layout.ti_i, layout.K, layout.K);
        auto          R      = input(u, series_.R());
        auto          L      = input(u, series_.L());
        auto          G      = input(u, shunt_.G());
        auto          C      = input(u, shunt_.C());

        auto XL = omega * L;
        auto BC = omega * C;

        auto TvRtrack = matrix(layout.K,
                               layout.K,
                               [&](IdxT i, IdxT m)
                               {
                                 V value = -amd[m] * TvR(i, m) + bmd[m] * TvI(i, m);
                                 for (IdxT k = 0; k < layout.K; ++k)
                                 {
                                   const auto delta  = (i == k) ? 1.0 : 0.0;
                                   const auto Ar     = Alpha(i, k) - am[m] * delta;
                                   const auto Br     = Beta(i, k) - bm[m] * delta;
                                   value            += Alphad(i, k) * TvR(k, m) - Betad(i, k) * TvI(k, m)
                                            + Ar * TvRd(k, m) - Br * TvId(k, m);
                                 }
                                 return value;
                               });
        auto TvItrack = matrix(layout.K,
                               layout.K,
                               [&](IdxT i, IdxT m)
                               {
                                 V value = -amd[m] * TvI(i, m) - bmd[m] * TvR(i, m);
                                 for (IdxT k = 0; k < layout.K; ++k)
                                 {
                                   const auto delta  = (i == k) ? 1.0 : 0.0;
                                   const auto Ar     = Alpha(i, k) - am[m] * delta;
                                   const auto Br     = Beta(i, k) - bm[m] * delta;
                                   value            += Alphad(i, k) * TvI(k, m) + Betad(i, k) * TvR(k, m)
                                            + Br * TvRd(k, m) + Ar * TvId(k, m);
                                 }
                                 return value;
                               });
        auto TiRtrack = matrix(layout.K,
                               layout.K,
                               [&](IdxT i, IdxT m)
                               {
                                 V value = -amd[m] * TiR(i, m) - bmd[m] * TiI(i, m);
                                 for (IdxT k = 0; k < layout.K; ++k)
                                 {
                                   const auto delta  = (i == k) ? 1.0 : 0.0;
                                   const auto Al     = Alpha(k, i) - am[m] * delta;
                                   const auto Bl     = Beta(k, i) - bm[m] * delta;
                                   value            += Alphad(k, i) * TiR(k, m) + Betad(k, i) * TiI(k, m)
                                            + Al * TiRd(k, m) + Bl * TiId(k, m);
                                 }
                                 return value;
                               });
        auto TiItrack = matrix(layout.K,
                               layout.K,
                               [&](IdxT i, IdxT m)
                               {
                                 V value = -amd[m] * TiI(i, m) + bmd[m] * TiR(i, m);
                                 for (IdxT k = 0; k < layout.K; ++k)
                                 {
                                   const auto delta  = (i == k) ? 1.0 : 0.0;
                                   const auto Al     = Alpha(k, i) - am[m] * delta;
                                   const auto Bl     = Beta(k, i) - bm[m] * delta;
                                   value            += Alphad(k, i) * TiI(k, m) - Betad(k, i) * TiR(k, m)
                                            - Bl * TiRd(k, m) + Al * TiId(k, m);
                                 }
                                 return value;
                               });
        auto aNorm    = vector(layout.K,
                            [&](IdxT m)
                            {
                              V value = -1.0;
                              for (IdxT i = 0; i < layout.K; ++i)
                              {
                                value += TiR(i, m) * TvR(i, m) + TiI(i, m) * TvI(i, m);
                              }
                              return value;
                            });
        auto bNorm    = vector(layout.K,
                            [&](IdxT m)
                            {
                              V value = 0.0;
                              for (IdxT i = 0; i < layout.K; ++i)
                              {
                                value += TiR(i, m) * TvI(i, m) - TiI(i, m) * TvR(i, m);
                              }
                              return value;
                            });

        equation(Ares)   = Alpha * Alpha - Beta * Beta - R * G + XL * BC;
        equation(Bres)   = Alpha * Beta + Beta * Alpha - R * BC - XL * G;
        equation(TvRres) = TvRtrack;
        equation(TvIres) = TvItrack;
        equation(TiRres) = TiRtrack;
        equation(TiIres) = TiItrack;
        equation(ares)   = aNorm;
        equation(bres)   = bNorm;
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
        return this->template evaluateInputDifferentialJacobian<Gamma<ScalarT, IdxT>>();
      }

      template class Gamma<double, long int>;
      template class Gamma<double, size_t>;
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
