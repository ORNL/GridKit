/**
 * @file SeriesReduction.cpp
 *
 * @brief Initialization, equation residual, and sparse-jet Jacobian for the
 * series reduction.
 *
 */

#include "SeriesReduction.hpp"

#include <complex>
#include <stdexcept>

#include <Eigen/LU>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename ScalarT, typename IdxT>
      SeriesReduction<ScalarT, IdxT>::SeriesReduction(
          const SeriesView<ScalarT, IdxT>& full, MapT map)
        : full_(full),
          map_(std::move(map)),
          layout_(map_.conductors, map_.terminals)
      {
      }

      template <typename ScalarT, typename IdxT>
      int SeriesReduction<ScalarT, IdxT>::initialize()
      {
        using ComplexT = std::complex<ScalarT>;
        using MatrixXc = Eigen::Matrix<ComplexT, Eigen::Dynamic, Eigen::Dynamic>;

        auto eigenIndex = [](IdxT i)
        { return static_cast<Eigen::Index>(i); };

        const Layout& layout = layout_;
        const ScalarT omega  = static_cast<ScalarT>(coordinate());
        if (!(omega > ScalarT{0}))
        {
          throw std::runtime_error(
              "Series reduction initialization requires a positive frequency");
        }

        auto R    = input(full_.R());
        auto L    = input(full_.L());
        auto Drs  = state(layout.dr, layout.C, layout.P);
        auto Dis  = state(layout.di, layout.C, layout.P);
        auto Rred = state(layout.rr, layout.P, layout.P);
        auto Lred = state(layout.lr, layout.P, layout.P);

        // The augmented block solve [Z, -E; E^T, 0][D; Zred] = [0; I]
        // states the elimination directly, so the residual vanishes at
        // this point.
        const IdxT n      = layout.C + layout.P;
        MatrixXc   system = MatrixXc::Zero(eigenIndex(n), eigenIndex(n));
        MatrixXc   rhs    = MatrixXc::Zero(eigenIndex(n), eigenIndex(layout.P));

        for (IdxT i = 0; i < layout.C; ++i)
        {
          for (IdxT j = 0; j < layout.C; ++j)
          {
            system(eigenIndex(i), eigenIndex(j)) =
                ComplexT{R(i, j), omega * L(i, j)};
          }
        }
        for (IdxT c = 0; c < layout.C; ++c)
        {
          const IdxT t = map_.terminal[static_cast<size_t>(c)];
          if (t != MapT::kGrounded)
          {
            system(eigenIndex(c), eigenIndex(layout.C + t)) = ComplexT{-1.0, 0.0};
            system(eigenIndex(layout.C + t), eigenIndex(c)) = ComplexT{1.0, 0.0};
          }
        }
        for (IdxT p = 0; p < layout.P; ++p)
        {
          rhs(eigenIndex(layout.C + p), eigenIndex(p)) = ComplexT{1.0, 0.0};
        }

        Eigen::FullPivLU<MatrixXc> lu(system);
        if (!lu.isInvertible())
        {
          throw std::runtime_error(
              "Series reduction requires an invertible constrained "
              "impedance block");
        }
        const MatrixXc solution = lu.solve(rhs);

        for (IdxT c = 0; c < layout.C; ++c)
        {
          for (IdxT p = 0; p < layout.P; ++p)
          {
            const auto value = solution(eigenIndex(c), eigenIndex(p));
            Drs(c, p)        = value.real();
            Dis(c, p)        = value.imag();
          }
        }
        for (IdxT t = 0; t < layout.P; ++t)
        {
          for (IdxT p = 0; p < layout.P; ++p)
          {
            const auto value = solution(eigenIndex(layout.C + t), eigenIndex(p));
            Rred(t, p)       = value.real();
            Lred(t, p)       = value.imag() / omega;
          }
        }

        return 0;
      }

      template <typename ScalarT, typename IdxT>
      template <typename V>
      __attribute__((always_inline)) inline void
      SeriesReduction<ScalarT, IdxT>::evaluateResidualKernel(V* y,
                                                             V*,
                                                             V* u,
                                                             V* f) const
      {
        const Layout& layout = layout_;
        const auto    omega  = static_cast<ScalarT>(coordinate());
        auto          Drs    = slice(y, layout.dr, layout.C, layout.P);
        auto          Dis    = slice(y, layout.di, layout.C, layout.P);
        auto          Rred   = slice(y, layout.rr, layout.P, layout.P);
        auto          Lred   = slice(y, layout.lr, layout.P, layout.P);
        auto          FDr    = slice(f, layout.dr, layout.C, layout.P);
        auto          FDi    = slice(f, layout.di, layout.C, layout.P);
        auto          FRr    = slice(f, layout.rr, layout.P, layout.P);
        auto          FRi    = slice(f, layout.lr, layout.P, layout.P);
        auto          R      = input(u, full_.R());
        auto          L      = input(u, full_.L());

        // Z D - E Zred = 0 with Z = R + j omega L, D = Dr + j Di, and
        // Zred = Rred + j omega Lred. E has one unit entry per bundled
        // conductor, so its product is an index lookup, not arithmetic.
        for (IdxT c = 0; c < layout.C; ++c)
        {
          const IdxT t = map_.terminal[static_cast<size_t>(c)];
          for (IdxT p = 0; p < layout.P; ++p)
          {
            V real{};
            V imag{};
            for (IdxT j = 0; j < layout.C; ++j)
            {
              real += R(c, j) * Drs(j, p) - omega * L(c, j) * Dis(j, p);
              imag += R(c, j) * Dis(j, p) + omega * L(c, j) * Drs(j, p);
            }
            if (t != MapT::kGrounded)
            {
              real -= Rred(t, p);
              imag -= omega * Lred(t, p);
            }
            FDr(c, p) = real;
            FDi(c, p) = imag;
          }
        }

        // E^T D - I = 0: terminal currents are bundle-member sums.
        for (IdxT t = 0; t < layout.P; ++t)
        {
          for (IdxT p = 0; p < layout.P; ++p)
          {
            V real{};
            V imag{};
            for (const IdxT c : map_.members[static_cast<size_t>(t)])
            {
              real += Drs(c, p);
              imag += Dis(c, p);
            }
            FRr(t, p) = real - static_cast<ScalarT>(t == p ? 1 : 0);
            FRi(t, p) = imag;
          }
        }
      }

      template <typename ScalarT, typename IdxT>
      void SeriesReduction<ScalarT, IdxT>::evaluateInternalResidual(
          ScalarT* y, ScalarT* yp, ScalarT* u, ScalarT* f) const
      {
        evaluateResidualKernel<ScalarT>(y, yp, u, f);
      }

      template <typename ScalarT, typename IdxT>
      int SeriesReduction<ScalarT, IdxT>::evaluateResidual()
      {
        evaluateInternalResidual(y_, yp_, nullptr, f_);
        return 0;
      }

      template <typename ScalarT, typename IdxT>
      int SeriesReduction<ScalarT, IdxT>::evaluateJacobian()
      {
        return this->template evaluateInputAlgebraicJacobian<
            SeriesReduction<ScalarT, IdxT>>();
      }

      template class SeriesReduction<double, long int>;
      template class SeriesReduction<double, size_t>;
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
