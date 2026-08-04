/**
 * @file VectorFitting.hpp
 *
 * @brief Rational approximation of sampled frequency responses posed as
 *        constrained optimization.
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#pragma once

#include <complex>
#include <string>
#include <vector>

#include <GridKit/ScalarTraits.hpp>
#include <GridKit/Solver/Optimization/Rational/RationalModel.hpp>
#include <GridKit/Solver/Optimization/Rational/SampledResponse.hpp>

namespace GridKit
{
  namespace Optimization
  {
    /// Optional terms alongside the pole-residue sum.
    enum class RationalTerms
    {
      NONE,     ///< Strictly proper: H(s) = sum of R_q / (s - p_q).
      CONSTANT, ///< Adds D.
      LINEAR    ///< Adds D + sE.
    };

    /// Named weighting schemes; see README for the conventions.
    enum class Weighting
    {
      UNIFORM,                   ///< Absolute error.
      INVERSE_MAGNITUDE,         ///< Balanced emphasis.
      INVERSE_SQUARED_MAGNITUDE, ///< Relative error.
      CUSTOM                     ///< Caller-supplied per-sample weights.
    };

    /// Initial pole placement over the sampled band.
    enum class PoleSeeding
    {
      LOG_COMPLEX, ///< Weakly damped conjugate pairs, log-spaced.
      LOG_REAL,    ///< Real poles, log-spaced.
      USER         ///< Caller-supplied starting poles.
    };

    /**
     * @brief Relaxed vector fitting in the real-augmented formulation.
     *
     * Each pole-relocation pass and the coefficient identification are
     * convex least squares problems solved in closed form by orthogonal
     * block elimination; the optional refinement stage is a
     * bound-constrained nonlinear program solved through the
     * Optimization family. See README.md for the algorithm
     * specification.
     */
    template <typename scalar_type, typename index_type>
    class VectorFitting
    {
    public:
      using ScalarT  = scalar_type;
      using IdxT     = index_type;
      using RealT    = typename GridKit::ScalarTraits<ScalarT>::RealT;
      using ComplexT = std::complex<RealT>;

      /// Algorithm options. Defaults reproduce relaxed vector fitting.
      struct Parameters
      {
        IdxT          pole_count = 10;
        RationalTerms terms      = RationalTerms::CONSTANT;
        Weighting     weighting  = Weighting::UNIFORM;

        /// Per-sample weights. Required with Weighting::CUSTOM and
        /// rejected otherwise; any mismatch is a hard error.
        std::vector<RealT> weights;

        /// Starting poles. Required with PoleSeeding::USER and rejected
        /// otherwise; any mismatch is a hard error.
        std::vector<ComplexT> starting_poles;

        PoleSeeding seeding        = PoleSeeding::LOG_COMPLEX;
        IdxT        max_iterations = 30;

        /// Largest nearest-match relative pole movement that counts a
        /// relocation pass as converged. The default sits above the
        /// noise floor of the equilibrated least squares solves.
        RealT pole_shift_tolerance = 1.0e-6;
        RealT stability_margin     = 0.0;

        /// Lowest-order search over the pole count.
        struct OrderSearch
        {
          bool  enabled        = false;
          IdxT  min_poles      = 2;
          IdxT  max_poles      = 30;
          RealT target_rel_rms = 1.0e-3;
        } order_search;

        /// Deterministic perturbation restarts; zero restarts disables.
        struct Restarts
        {
          IdxT  max_restarts    = 0;
          RealT trigger_rel_rms = 0.05;
        } restarts;

        /// Stage C nonlinear refinement toggle.
        bool refine = false;
      };

      /// Record of one pole-relocation pass.
      struct IterationInfo
      {
        IdxT  iteration{0};
        RealT max_pole_shift{0.0};
      };

      /// Statistics since the last fit() call.
      struct Stats
      {
        std::vector<IterationInfo> iterations;
        std::vector<RealT>         channel_rel_rms;
        RealT                      final_rel_rms{0.0};
        RealT                      final_abs_rms{0.0};
        bool                       converged{false};
        IdxT                       restarts_used{0};
        IdxT                       order_selected{0};

        std::string report() const;
      };

      /**
       * @brief Bind the fitter to a sampled response.
       *
       * @pre The samples outlive the fitter and are not modified while a
       *      fit is in progress
       */
      explicit VectorFitting(const SampledResponse<ScalarT, IdxT>& samples);

      /**
       * @brief Fit a rational approximation to the sampled response.
       *
       * @param[out] model  Fitted poles, residues, and polynomial terms
       * @param[in]  params Algorithm options
       * @return 0 on success, 1 when the relocation iteration limit was
       *         reached (result usable), 2 when the order search ended
       *         without meeting its target (best model returned),
       *         negative on a hard failure (inconsistent options,
       *         dimension mismatch, optimizer or eigensolver error)
       *
       * @pre The sample set has strictly increasing, positive frequencies
       * @post Every pole satisfies the stability margin
       */
      [[nodiscard("May fail. Check error code.")]]
      int fit(RationalModel<ScalarT, IdxT>& model, Parameters params = {});

      const Stats& getStats() const
      {
        return stats_;
      }

    private:
      const SampledResponse<ScalarT, IdxT>& samples_;
      Stats                                 stats_;
    };
  } // namespace Optimization
} // namespace GridKit
