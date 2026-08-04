/**
 * @file ModalDecomposition.hpp
 *
 * @brief Identity-preserving modal observation of the propagation
 *        matrix along a frequency sweep.
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#pragma once

#include <complex>
#include <vector>

#include <GridKit/Model/EMT/Signal.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      /**
       * @brief Spectral decomposition of Gamma with mode identity owned
       *        across samples.
       *
       * The frequency sweep integrates only quantities analytic in
       * omega; individual eigenpairs are not among them, so they are
       * observed here, outside the DAE. Each decompose() call factors
       * the current Gamma exactly and assigns mode identity by
       * biorthogonal overlap with the previous sample, so identity
       * follows the eigenvectors and survives eigenvalue crossings.
       *
       * Modes whose eigenvalues fall within the relative cluster gap
       * form a cluster: only their invariant subspace is defined, so
       * the reported basis is the orthonormal span of the cluster
       * rotated by the polar factor of its overlap with the previous
       * frame, the closest admissible continuation. A simple mode is
       * the one-by-one case of the same rule, its polar factor a
       * phase, so a single alignment law covers every mode.
       *
       * Ti^H Tv = I holds to machine precision at every sample because
       * the duals are computed, never transported. Presenting a
       * cluster through its own eigenvalues bounds the defect of
       * Tv diag(e^{-lambda l}) Ti^H against e^{-Gamma l} by
       * O(l |gamma| cluster_gap); at exact degeneracy it vanishes.
       *
       * Results live in one flat buffer described by Signal views, so
       * a variable monitor can observe them exactly as it observes
       * element states.
       */
      template <typename scalar_type, typename index_type>
      class ModalDecomposition
      {
      public:
        using ScalarT  = scalar_type;
        using IdxT     = index_type;
        using RealT    = typename GridKit::ScalarTraits<ScalarT>::RealT;
        using ComplexT = std::complex<RealT>;
        using SignalT  = Signal<IdxT>;

        struct Parameters
        {
          /// Relative eigenvalue distance below which modes form one
          /// cluster and only their common subspace is meaningful.
          RealT cluster_gap = 1.0e-8;
        };

        /**
         * @brief Bind the mode count and the line length the modal
         *        responses are evaluated for.
         *
         * @pre conductor_count >= 1 and line_length > 0; the
         *      constructor throws std::runtime_error otherwise
         */
        ModalDecomposition(IdxT       conductor_count,
                           RealT      line_length,
                           Parameters parameters = {});

        /**
         * @brief Observe Gamma = Alpha + j Beta at one sample.
         *
         * The first call fixes the canonical mode order and gauge;
         * every later call preserves identity and continuity against
         * the previous sample. Samples must arrive in sweep order.
         *
         * @param[in] omega  Angular frequency of the sample, positive
         * @param[in] alpha  Re(Gamma), row-major K by K
         * @param[in] beta   Im(Gamma), row-major K by K
         * @return 0 on success, -1 on a nonpositive frequency or a
         *         failed eigensolve
         */
        [[nodiscard("May fail. Check error code.")]]
        int decompose(RealT omega, const ScalarT* alpha, const ScalarT* beta);

        /// Forget the tracked frame so the next decompose() starts a
        /// fresh sweep with the canonical order and gauge.
        void reset();

        IdxT modeCount() const
        {
          return layout_.K;
        }

        /// Storage the signal views index into; valid after the first
        /// successful decompose().
        const RealT* data() const
        {
          return values_.data();
        }

        SignalT tvReal() const
        {
          return {layout_.tv_r, layout_.K, layout_.K};
        }

        SignalT tvImag() const
        {
          return {layout_.tv_i, layout_.K, layout_.K};
        }

        SignalT tiReal() const
        {
          return {layout_.ti_r, layout_.K, layout_.K};
        }

        SignalT tiImag() const
        {
          return {layout_.ti_i, layout_.K, layout_.K};
        }

        /// Modal attenuation constants Re(lambda_m).
        SignalT alphaM() const
        {
          return {layout_.a, layout_.K, 1};
        }

        /// Modal phase constants Im(lambda_m).
        SignalT betaM() const
        {
          return {layout_.b, layout_.K, 1};
        }

        /// Propagation factors e^{-lambda_m l}.
        SignalT hReal() const
        {
          return {layout_.h_r, layout_.K, 1};
        }

        SignalT hImag() const
        {
          return {layout_.h_i, layout_.K, 1};
        }

        /// Modal phase delays l Im(lambda_m) / omega.
        SignalT tau() const
        {
          return {layout_.tau, layout_.K, 1};
        }

      private:
        struct Layout
        {
          explicit Layout(IdxT conductor_count)
            : K(conductor_count),
              KK(K * K),
              tv_r(0),
              tv_i(KK),
              ti_r(2 * KK),
              ti_i(3 * KK),
              a(4 * KK),
              b(a + K),
              h_r(b + K),
              h_i(h_r + K),
              tau(h_i + K),
              n(tau + K)
          {
          }

          const IdxT K;
          const IdxT KK;
          const IdxT tv_r;
          const IdxT tv_i;
          const IdxT ti_r;
          const IdxT ti_i;
          const IdxT a;
          const IdxT b;
          const IdxT h_r;
          const IdxT h_i;
          const IdxT tau;
          const IdxT n;
        };

        const Parameters parameters_;
        const Layout     layout_;
        const RealT      length_;

        std::vector<RealT> values_;

        /// Previous frame, the reference every alignment is made
        /// against. Empty until the first decompose().
        std::vector<ComplexT> tv_previous_;
        std::vector<ComplexT> ti_previous_;
        std::vector<ComplexT> lambda_previous_;
      };
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
