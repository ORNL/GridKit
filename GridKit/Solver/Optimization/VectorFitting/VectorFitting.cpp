/**
 * @file VectorFitting.cpp
 *
 * @brief Relaxed vector fitting in the real-augmented formulation.
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#include "VectorFitting.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <utility>

#include <GridKit/Solver/Optimization/ConstrainedLeastSquares.hpp>
#include <GridKit/Solver/Optimization/IpoptSolver.hpp>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace
{
  constexpr double SEED_DAMPING_RATIO    = 0.01;
  constexpr double SIGMA_SCALE_FLOOR     = 1.0e-30;
  constexpr double SHIFT_NORM_FLOOR      = 1.0e-15;
  constexpr double REAL_EIGENVALUE_TOL   = 1.0e-12;
  constexpr double WEIGHT_FLOOR_RATIO    = 1.0e-8;
  constexpr double COLUMN_NORM_FLOOR     = 1.0e-30;
  constexpr double DEGENERATE_RELOCATION = 1.0e-14;
  constexpr double QP_TOLERANCE          = 1.0e-10;
  constexpr size_t QP_MAX_ITERATIONS     = 500;
  constexpr double REFINE_TOLERANCE      = 1.0e-8;
  constexpr size_t REFINE_MAX_ITERATIONS = 300;
  constexpr double RESTART_GOLDEN_FRACT  = 0.6180339887;

  using GridKit::Optimization::ConstrainedLeastSquares;
  using GridKit::Optimization::HessianMode;
  using GridKit::Optimization::IpoptSolver;
  using GridKit::Optimization::OptimizationModel;
  using GridKit::Optimization::RationalModel;
  using GridKit::Optimization::RationalTerms;
  using GridKit::Optimization::SampledResponse;
  using GridKit::Optimization::Weighting;

  size_t termColumns(RationalTerms terms)
  {
    if (terms == RationalTerms::NONE)
    {
      return 0;
    }
    return terms == RationalTerms::CONSTANT ? 1 : 2;
  }

  /**
   * @brief Real-coefficient basis values at s = j omega for a canonical
   *        pole set: a real pole q contributes 1/(s - p_q); a conjugate
   *        pair occupying slots q and q+1 contributes
   *        1/(s - p) + 1/(s - conj(p)) and j/(s - p) - j/(s - conj(p)).
   */
  template <typename real_type>
  void fillBasis(const std::vector<std::complex<real_type>>& poles,
                 real_type                                   omega,
                 std::vector<std::complex<real_type>>&       phi)
  {
    using ComplexT = std::complex<real_type>;

    const ComplexT s{real_type{0}, omega};
    const ComplexT unit{real_type{0}, real_type{1}};

    size_t q = 0;
    while (q < poles.size())
    {
      if (poles[q].imag() == real_type{0})
      {
        phi[q]  = real_type{1} / (s - poles[q]);
        q      += 1;
        continue;
      }
      const ComplexT first   = real_type{1} / (s - poles[q]);
      const ComplexT second  = real_type{1} / (s - std::conj(poles[q]));
      phi[q]                 = first + second;
      phi[q + 1]             = unit * (first - second);
      q                     += 2;
    }
  }

  /**
   * @brief Flip eigenvalues into the stable half plane and order them
   *        canonically: real poles ascending, then conjugate pairs by
   *        ascending imaginary part with the positive member first.
   */
  template <typename real_type>
  std::vector<std::complex<real_type>>
  canonicalizePoles(const std::vector<std::complex<real_type>>& raw,
                    real_type                                   margin)
  {
    using ComplexT = std::complex<real_type>;

    std::vector<real_type> reals;
    std::vector<ComplexT>  pairs;

    for (const auto& value : raw)
    {
      const real_type a =
          std::min(-std::abs(value.real()), -margin);
      const real_type scale =
          real_type{REAL_EIGENVALUE_TOL} * (real_type{1} + std::abs(value));

      if (std::abs(value.imag()) <= scale)
      {
        reals.push_back(a);
      }
      else if (value.imag() > real_type{0})
      {
        pairs.push_back(ComplexT{a, value.imag()});
      }
    }

    std::sort(reals.begin(), reals.end());
    std::sort(pairs.begin(), pairs.end(), [](ComplexT lhs, ComplexT rhs)
              { return lhs.imag() < rhs.imag(); });

    std::vector<ComplexT> poles;
    poles.reserve(reals.size() + 2 * pairs.size());
    for (const auto a : reals)
    {
      poles.push_back(ComplexT{a, real_type{0}});
    }
    for (const auto& pole : pairs)
    {
      poles.push_back(pole);
      poles.push_back(std::conj(pole));
    }
    return poles;
  }

  /**
   * @brief Largest relative distance from a relocated pole to its
   *        nearest previous pole. Nearest-match distances are immune to
   *        the re-sorting canonicalizePoles performs, which a positional
   *        comparison would misread as a large shift.
   */
  template <typename real_type>
  real_type
  maxRelativeShift(const std::vector<std::complex<real_type>>& previous,
                   const std::vector<std::complex<real_type>>& current)
  {
    real_type shift = real_type{0};
    for (const auto& pole : current)
    {
      real_type nearest = std::numeric_limits<real_type>::infinity();
      for (const auto& reference : previous)
      {
        const real_type scale =
            std::abs(reference) + real_type{SHIFT_NORM_FLOOR};
        nearest = std::min(nearest, std::abs(pole - reference) / scale);
      }
      if (!previous.empty())
      {
        shift = std::max(shift, nearest);
      }
    }
    return shift;
  }
} // namespace

namespace GridKit
{
  namespace Optimization
  {
    namespace Detail
    {
      /**
       * @brief Solve one convex quadratic program through the family
       *        backend with the exact precomputed Hessian.
       *
       * The columns are equilibrated to unit weighted norm before the
       * solve: the basis magnitudes span the sampled band's decades, and
       * the precomputed normal-equations Hessian squares that spread
       * into its condition number. The solution is mapped back to the
       * caller's variables afterwards.
       *
       * @return 0 when the solver certifies optimal or acceptable
       *         convergence; -2 otherwise, including an iteration-limit
       *         stop, which carries no convergence certificate.
       */
      template <typename scalar_type, typename index_type>
      int solveQuadraticProgram(
          typename ConstrainedLeastSquares<scalar_type,
                                           index_type>::Problem problem,
          std::vector<typename ConstrainedLeastSquares<
              scalar_type,
              index_type>::RealT>&                              x)
      {
        using SolverT = IpoptSolver<scalar_type, index_type>;
        using RealT   = typename ConstrainedLeastSquares<scalar_type,
                                                         index_type>::RealT;

        const auto variables  = static_cast<size_t>(problem.variable_count);
        const auto residuals  = static_cast<size_t>(problem.residual_count);
        const auto equalities = static_cast<size_t>(problem.equality_count);

        std::vector<RealT> scale(variables, RealT{1});
        for (size_t c = 0; c < variables; ++c)
        {
          RealT norm = RealT{0};
          for (size_t m = 0; m < residuals; ++m)
          {
            const RealT w =
                problem.weights.empty() ? RealT{1} : problem.weights[m];
            const RealT entry  = problem.a[m * variables + c];
            norm              += w * entry * entry;
          }
          norm = std::sqrt(norm);
          if (norm > RealT{COLUMN_NORM_FLOOR})
          {
            scale[c] = RealT{1} / norm;
          }
        }
        for (size_t m = 0; m < residuals; ++m)
        {
          for (size_t c = 0; c < variables; ++c)
          {
            problem.a[m * variables + c] *= scale[c];
          }
        }
        for (size_t i = 0; i < equalities; ++i)
        {
          for (size_t c = 0; c < variables; ++c)
          {
            problem.a_eq[i * variables + c] *= scale[c];
          }
        }

        try
        {
          ConstrainedLeastSquares<scalar_type, index_type> model(
              std::move(problem));

          typename SolverT::Parameters options;
          options.tolerance = QP_TOLERANCE;
          options.max_iterations =
              static_cast<index_type>(QP_MAX_ITERATIONS);
          options.hessian = HessianMode::EXACT;

          SolverT   solver(&model, options);
          const int status = solver.solve(x);
          if (status < 0 || status > 1)
          {
            return -2;
          }
          for (size_t c = 0; c < variables; ++c)
          {
            x[c] *= scale[c];
          }
          return 0;
        }
        catch (const std::exception&)
        {
          return -2;
        }
      }

      /**
       * @brief One relaxed relocation pass: solve the sigma quadratic
       *        program with the relaxation equality, then relocate the
       *        poles to the zeros of sigma via the real companion
       *        eigenproblem.
       *
       * A solution with numerator coefficients vanishing against the
       * relaxation constant makes sigma effectively constant: the
       * companion update is below machine noise and the poles do not
       * move for structural reasons. Such a pass reports @p degenerate
       * and must not count as a convergence certificate.
       */
      template <typename scalar_type, typename index_type>
      int relocatePoles(
          const SampledResponse<scalar_type, index_type>&           samples,
          const std::vector<typename SampledResponse<
              scalar_type,
              index_type>::RealT>&                                  sample_weights,
          RationalTerms                                             terms,
          typename SampledResponse<scalar_type, index_type>::RealT  margin,
          std::vector<typename SampledResponse<
              scalar_type,
              index_type>::ComplexT>&                               poles,
          typename SampledResponse<scalar_type, index_type>::RealT& shift,
          bool&                                                     degenerate)
      {
        using RealT    = typename SampledResponse<scalar_type,
                                                  index_type>::RealT;
        using ComplexT = typename SampledResponse<scalar_type,
                                                  index_type>::ComplexT;
        using ProblemT = typename ConstrainedLeastSquares<
            scalar_type,
            index_type>::Problem;

        const auto channels =
            static_cast<size_t>(samples.rows) * static_cast<size_t>(samples.cols);
        const auto sample_count = samples.omega.size();
        const auto order        = poles.size();
        const auto numerator    = order + termColumns(terms);
        const auto variables    = channels * numerator + order + 1;
        const auto rows         = 2 * sample_count * channels;

        ProblemT problem;
        problem.variable_count = static_cast<index_type>(variables);
        problem.residual_count = static_cast<index_type>(rows);
        problem.equality_count = 1;
        problem.a.assign(rows * variables, RealT{0});
        problem.b.assign(rows, RealT{0});
        problem.weights.assign(rows, RealT{1});
        problem.a_eq.assign(variables, RealT{0});
        problem.b_eq.assign(1, static_cast<RealT>(sample_count));

        const size_t          sigma_offset = channels * numerator;
        std::vector<ComplexT> phi(order);

        for (size_t m = 0; m < sample_count; ++m)
        {
          fillBasis(poles, samples.omega[m], phi);

          for (size_t k = 0; k < order; ++k)
          {
            problem.a_eq[sigma_offset + k] += phi[k].real();
          }
          problem.a_eq[sigma_offset + order] += RealT{1};

          for (size_t ch = 0; ch < channels; ++ch)
          {
            const ComplexT response = samples.response[m * channels + ch];
            const size_t   row_re   = 2 * (m * channels + ch);
            const size_t   row_im   = row_re + 1;
            RealT*         re       = &problem.a[row_re * variables];
            RealT*         im       = &problem.a[row_im * variables];

            for (size_t k = 0; k < order; ++k)
            {
              re[ch * numerator + k] = phi[k].real();
              im[ch * numerator + k] = phi[k].imag();

              const ComplexT scaled = response * phi[k];
              re[sigma_offset + k]  = -scaled.real();
              im[sigma_offset + k]  = -scaled.imag();
            }
            if (termColumns(terms) > 0)
            {
              re[ch * numerator + order] = RealT{1};
            }
            if (termColumns(terms) > 1)
            {
              im[ch * numerator + order + 1] = samples.omega[m];
            }
            re[sigma_offset + order] = -response.real();
            im[sigma_offset + order] = -response.imag();

            problem.weights[row_re] = sample_weights[m];
            problem.weights[row_im] = sample_weights[m];
          }
        }

        std::vector<RealT> solution;
        if (solveQuadraticProgram<scalar_type, index_type>(
                std::move(problem), solution)
            != 0)
        {
          return -2;
        }

        RealT sigma_scale = solution[sigma_offset + order];
        if (std::abs(sigma_scale) < RealT{SIGMA_SCALE_FLOOR})
        {
          sigma_scale = sigma_scale < RealT{0} ? -RealT{SIGMA_SCALE_FLOOR}
                                               : RealT{SIGMA_SCALE_FLOOR};
        }

        RealT gamma_norm = RealT{0};
        RealT pole_scale = RealT{0};
        for (size_t k = 0; k < order; ++k)
        {
          gamma_norm += solution[sigma_offset + k] * solution[sigma_offset + k];
          pole_scale  = std::max(pole_scale, std::abs(poles[k]));
        }
        gamma_norm = std::sqrt(gamma_norm);

        degenerate = RealT{2} * gamma_norm
                     < RealT{DEGENERATE_RELOCATION} * std::abs(sigma_scale)
                           * (RealT{1} + pole_scale);
        if (degenerate)
        {
          shift = RealT{0};
          return 0;
        }

        // Zeros of sigma: eigenvalues of A - b gamma^T / sigma_scale with
        // the pole matrix in its real block form. With b = [2, 0]^T per
        // pair the basis coefficients enter gamma unchanged: the block
        // realizes c' phi + c'' psi exactly.
        Eigen::MatrixXd companion =
            Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(order),
                                  static_cast<Eigen::Index>(order));
        Eigen::VectorXd b_vec =
            Eigen::VectorXd::Zero(static_cast<Eigen::Index>(order));
        Eigen::VectorXd gamma =
            Eigen::VectorXd::Zero(static_cast<Eigen::Index>(order));

        size_t q = 0;
        while (q < order)
        {
          const auto row = static_cast<Eigen::Index>(q);
          if (poles[q].imag() == RealT{0})
          {
            companion(row, row)  = poles[q].real();
            b_vec(row)           = 1.0;
            gamma(row)           = solution[sigma_offset + q];
            q                   += 1;
            continue;
          }
          companion(row, row)          = poles[q].real();
          companion(row, row + 1)      = poles[q].imag();
          companion(row + 1, row)      = -poles[q].imag();
          companion(row + 1, row + 1)  = poles[q].real();
          b_vec(row)                   = 2.0;
          gamma(row)                   = solution[sigma_offset + q];
          gamma(row + 1)               = solution[sigma_offset + q + 1];
          q                           += 2;
        }
        companion -= (b_vec * gamma.transpose()) / sigma_scale;

        const Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(companion);
        if (eigensolver.info() != Eigen::Success)
        {
          return -2;
        }

        std::vector<ComplexT> raw(order);
        for (size_t k = 0; k < order; ++k)
        {
          const auto value =
              eigensolver.eigenvalues()(static_cast<Eigen::Index>(k));
          raw[k] = ComplexT{value.real(), value.imag()};
        }

        auto relocated = canonicalizePoles(raw, margin);
        shift          = maxRelativeShift(poles, relocated);
        poles          = std::move(relocated);
        return 0;
      }

      /**
       * @brief Identify residues and polynomial terms with the poles
       *        frozen: one unconstrained convex quadratic program per
       *        channel over the shared real basis.
       */
      template <typename scalar_type, typename index_type>
      int identifyCoefficients(
          const SampledResponse<scalar_type, index_type>& samples,
          const std::vector<typename SampledResponse<
              scalar_type,
              index_type>::RealT>&                        sample_weights,
          const std::vector<typename SampledResponse<
              scalar_type,
              index_type>::ComplexT>&                     poles,
          RationalTerms                                   terms,
          RationalModel<scalar_type, index_type>&         model)
      {
        using RealT    = typename SampledResponse<scalar_type,
                                                  index_type>::RealT;
        using ComplexT = typename SampledResponse<scalar_type,
                                                  index_type>::ComplexT;
        using ProblemT = typename ConstrainedLeastSquares<
            scalar_type,
            index_type>::Problem;

        const auto channels =
            static_cast<size_t>(samples.rows) * static_cast<size_t>(samples.cols);
        const auto sample_count = samples.omega.size();
        const auto order        = poles.size();
        const auto numerator    = order + termColumns(terms);

        model.rows  = samples.rows;
        model.cols  = samples.cols;
        model.poles = poles;
        model.residues.assign(order * channels, ComplexT{});
        model.d.clear();
        model.e.clear();
        if (termColumns(terms) > 0)
        {
          model.d.assign(channels, RealT{0});
        }
        if (termColumns(terms) > 1)
        {
          model.e.assign(channels, RealT{0});
        }

        std::vector<ComplexT>              phi(order);
        std::vector<std::vector<ComplexT>> basis(sample_count);
        for (size_t m = 0; m < sample_count; ++m)
        {
          fillBasis(poles, samples.omega[m], phi);
          basis[m] = phi;
        }

        for (size_t ch = 0; ch < channels; ++ch)
        {
          ProblemT problem;
          problem.variable_count = static_cast<index_type>(numerator);
          problem.residual_count =
              static_cast<index_type>(2 * sample_count);
          problem.equality_count = 0;
          problem.a.assign(2 * sample_count * numerator, RealT{0});
          problem.b.assign(2 * sample_count, RealT{0});
          problem.weights.assign(2 * sample_count, RealT{1});

          for (size_t m = 0; m < sample_count; ++m)
          {
            const ComplexT response = samples.response[m * channels + ch];
            RealT*         re       = &problem.a[(2 * m) * numerator];
            RealT*         im       = &problem.a[(2 * m + 1) * numerator];

            for (size_t k = 0; k < order; ++k)
            {
              re[k] = basis[m][k].real();
              im[k] = basis[m][k].imag();
            }
            if (termColumns(terms) > 0)
            {
              re[order] = RealT{1};
            }
            if (termColumns(terms) > 1)
            {
              im[order + 1] = samples.omega[m];
            }
            problem.b[2 * m]           = response.real();
            problem.b[2 * m + 1]       = response.imag();
            problem.weights[2 * m]     = sample_weights[m];
            problem.weights[2 * m + 1] = sample_weights[m];
          }

          std::vector<RealT> solution;
          if (solveQuadraticProgram<scalar_type, index_type>(
                  std::move(problem), solution)
              != 0)
          {
            return -2;
          }

          size_t q = 0;
          while (q < order)
          {
            if (poles[q].imag() == RealT{0})
            {
              model.residues[q * channels + ch] =
                  ComplexT{solution[q], RealT{0}};
              q += 1;
              continue;
            }
            const ComplexT residue{solution[q], solution[q + 1]};
            model.residues[q * channels + ch]        = residue;
            model.residues[(q + 1) * channels + ch]  = std::conj(residue);
            q                                       += 2;
          }
          if (termColumns(terms) > 0)
          {
            model.d[ch] = solution[order];
          }
          if (termColumns(terms) > 1)
          {
            model.e[ch] = solution[order + 1];
          }
        }
        return 0;
      }

      /**
       * @brief Unweighted error statistics over all samples and channels.
       */
      template <typename scalar_type, typename index_type>
      void computeStats(
          const SampledResponse<scalar_type, index_type>&         samples,
          const RationalModel<scalar_type, index_type>&           model,
          typename VectorFitting<scalar_type, index_type>::Stats& stats)
      {
        using RealT    = typename SampledResponse<scalar_type,
                                                  index_type>::RealT;
        using ComplexT = typename SampledResponse<scalar_type,
                                                  index_type>::ComplexT;

        const auto channels     = static_cast<size_t>(samples.rows) * static_cast<size_t>(samples.cols);
        const auto sample_count = samples.omega.size();

        RealT misfit    = RealT{0};
        RealT reference = RealT{0};
        stats.channel_rel_rms.assign(channels, RealT{0});
        std::vector<RealT> channel_reference(channels, RealT{0});

        for (size_t m = 0; m < sample_count; ++m)
        {
          for (size_t ch = 0; ch < channels; ++ch)
          {
            const auto row =
                static_cast<index_type>(ch) / samples.cols;
            const auto col =
                static_cast<index_type>(ch) % samples.cols;
            const ComplexT sample = samples.response[m * channels + ch];
            const ComplexT error =
                model.evaluate(samples.omega[m], row, col) - sample;
            misfit                    += std::norm(error);
            reference                 += std::norm(sample);
            stats.channel_rel_rms[ch] += std::norm(error);
            channel_reference[ch]     += std::norm(sample);
          }
        }

        const auto entries  = static_cast<RealT>(sample_count * channels);
        stats.final_abs_rms = std::sqrt(misfit / entries);
        stats.final_rel_rms =
            stats.final_abs_rms / std::max(std::sqrt(reference / entries), RealT{SHIFT_NORM_FLOOR});
        for (size_t ch = 0; ch < channels; ++ch)
        {
          stats.channel_rel_rms[ch] = std::sqrt(
              stats.channel_rel_rms[ch] / std::max(channel_reference[ch], RealT{SHIFT_NORM_FLOOR}));
        }
      }

      /**
       * @brief Stage C: bound-constrained refinement of the true weighted
       *        misfit over pole parameters and coefficients.
       *
       * Variables pack as [pole parameters, per-channel coefficients,
       * D, E]: a real pole carries its position a_q; a conjugate pair
       * carries (a_q, omega_q). Stability is the variable bound
       * a_q <= -margin and omega_q >= 0, so the program has no general
       * constraints and runs with a quasi-Newton Hessian.
       */
      template <typename scalar_type, typename index_type>
      class RefinementProblem final
        : public OptimizationModel<scalar_type, index_type>
      {
      public:
        using ScalarT  = scalar_type;
        using IdxT     = index_type;
        using RealT    = typename OptimizationModel<ScalarT, IdxT>::RealT;
        using ComplexT = std::complex<RealT>;

        RefinementProblem(
            const SampledResponse<ScalarT, IdxT>& samples,
            const std::vector<RealT>&             sample_weights,
            const RationalModel<ScalarT, IdxT>&   model,
            RationalTerms                         terms,
            RealT                                 margin)
          : samples_(samples),
            weights_(sample_weights),
            structure_(model.poles),
            terms_(terms),
            margin_(margin),
            channels_(static_cast<size_t>(model.rows) * static_cast<size_t>(model.cols)),
            order_(model.poles.size()),
            numerator_(model.poles.size() + termColumns(terms)),
            start_(model.poles.size() + (static_cast<size_t>(model.rows) * static_cast<size_t>(model.cols)) * (model.poles.size() + termColumns(terms)),
                   RealT{0})
        {
          size_t q = 0;
          while (q < order_)
          {
            start_[q] = structure_[q].real();
            if (structure_[q].imag() != RealT{0})
            {
              start_[q + 1]  = structure_[q].imag();
              q             += 2;
              continue;
            }
            q += 1;
          }
          for (size_t ch = 0; ch < channels_; ++ch)
          {
            size_t slot = 0;
            while (slot < order_)
            {
              const ComplexT residue =
                  model.residues[slot * channels_ + ch];
              start_[coefIndex(ch, slot)] = residue.real();
              if (structure_[slot].imag() != RealT{0})
              {
                start_[coefIndex(ch, slot + 1)]  = residue.imag();
                slot                            += 2;
                continue;
              }
              slot += 1;
            }
            if (termColumns(terms_) > 0)
            {
              start_[coefIndex(ch, order_)] = model.d[ch];
            }
            if (termColumns(terms_) > 1)
            {
              start_[coefIndex(ch, order_ + 1)] = model.e[ch];
            }
          }
        }

        IdxT variableCount() const override
        {
          return static_cast<IdxT>(start_.size());
        }

        IdxT constraintCount() const override
        {
          return 0;
        }

        int variableBounds(RealT* lower, RealT* upper) const override
        {
          const auto infinity = std::numeric_limits<RealT>::infinity();
          for (size_t k = 0; k < start_.size(); ++k)
          {
            lower[k] = -infinity;
            upper[k] = infinity;
          }
          size_t q = 0;
          while (q < order_)
          {
            upper[q] = -margin_;
            if (structure_[q].imag() != RealT{0})
            {
              lower[q + 1]  = RealT{0};
              q            += 2;
              continue;
            }
            q += 1;
          }
          return 0;
        }

        int constraintBounds(RealT*, RealT*) const override
        {
          return 0;
        }

        int startingPoint(RealT* x) const override
        {
          for (size_t k = 0; k < start_.size(); ++k)
          {
            x[k] = start_[k];
          }
          return 0;
        }

        int evaluateObjective(const RealT* x, RealT& value) const override
        {
          value                       = RealT{0};
          const auto            poles = buildPoles(x);
          std::vector<ComplexT> phi(order_);
          for (size_t m = 0; m < samples_.omega.size(); ++m)
          {
            fillBasis(poles, samples_.omega[m], phi);
            for (size_t ch = 0; ch < channels_; ++ch)
            {
              value += weights_[m] * std::norm(channelResidual(x, m, ch, phi));
            }
          }
          return 0;
        }

        int evaluateGradient(const RealT* x, RealT* gradient) const override
        {
          for (size_t k = 0; k < start_.size(); ++k)
          {
            gradient[k] = RealT{0};
          }

          const auto            poles = buildPoles(x);
          std::vector<ComplexT> phi(order_);
          std::vector<ComplexT> da(order_);
          std::vector<ComplexT> dw(order_);
          const ComplexT        unit{RealT{0}, RealT{1}};

          for (size_t m = 0; m < samples_.omega.size(); ++m)
          {
            const ComplexT s{RealT{0}, samples_.omega[m]};
            fillBasis(poles, samples_.omega[m], phi);

            size_t q = 0;
            while (q < order_)
            {
              if (structure_[q].imag() == RealT{0})
              {
                da[q]  = phi[q] * phi[q];
                q     += 1;
                continue;
              }
              const ComplexT first = RealT{1} / (s - poles[q]);
              const ComplexT second =
                  RealT{1} / (s - std::conj(poles[q]));
              const ComplexT sum = first * first + second * second;
              const ComplexT diff =
                  unit * (first * first - second * second);
              da[q]      = sum;
              da[q + 1]  = diff;
              dw[q]      = diff;
              dw[q + 1]  = -sum;
              q         += 2;
            }

            for (size_t ch = 0; ch < channels_; ++ch)
            {
              const ComplexT residual = channelResidual(x, m, ch, phi);
              const RealT    factor   = RealT{2} * weights_[m];
              const ComplexT weight   = factor * std::conj(residual);

              for (size_t k = 0; k < order_; ++k)
              {
                gradient[coefIndex(ch, k)] += (weight * phi[k]).real();
              }
              if (termColumns(terms_) > 0)
              {
                gradient[coefIndex(ch, order_)] += weight.real();
              }
              if (termColumns(terms_) > 1)
              {
                gradient[coefIndex(ch, order_ + 1)] +=
                    (weight * s).real();
              }

              size_t slot = 0;
              while (slot < order_)
              {
                const RealT alpha = x[coefIndex(ch, slot)];
                if (structure_[slot].imag() == RealT{0})
                {
                  gradient[slot] += (weight * (alpha * da[slot])).real();
                  slot           += 1;
                  continue;
                }
                const RealT beta = x[coefIndex(ch, slot + 1)];
                gradient[slot] +=
                    (weight * (alpha * da[slot] + beta * da[slot + 1]))
                        .real();
                gradient[slot + 1] +=
                    (weight * (alpha * dw[slot] + beta * dw[slot + 1]))
                        .real();
                slot += 2;
              }
            }
          }
          return 0;
        }

        int evaluateConstraints(const RealT*, RealT*) const override
        {
          return 0;
        }

        int evaluateJacobian(const RealT*, RealT*) const override
        {
          return 0;
        }

        /**
         * @brief Rebuild a rational model from a packed variable vector.
         */
        void unpack(const RealT*                  x,
                    RationalModel<ScalarT, IdxT>& model) const
        {
          size_t q = 0;
          while (q < order_)
          {
            if (structure_[q].imag() == RealT{0})
            {
              model.poles[q]  = ComplexT{x[q], RealT{0}};
              q              += 1;
              continue;
            }
            model.poles[q]      = ComplexT{x[q], x[q + 1]};
            model.poles[q + 1]  = std::conj(model.poles[q]);
            q                  += 2;
          }
          for (size_t ch = 0; ch < channels_; ++ch)
          {
            size_t slot = 0;
            while (slot < order_)
            {
              if (structure_[slot].imag() == RealT{0})
              {
                model.residues[slot * channels_ + ch] =
                    ComplexT{x[coefIndex(ch, slot)], RealT{0}};
                slot += 1;
                continue;
              }
              const ComplexT residue{x[coefIndex(ch, slot)],
                                     x[coefIndex(ch, slot + 1)]};
              model.residues[slot * channels_ + ch] = residue;
              model.residues[(slot + 1) * channels_ + ch] =
                  std::conj(residue);
              slot += 2;
            }
            if (termColumns(terms_) > 0)
            {
              model.d[ch] = x[coefIndex(ch, order_)];
            }
            if (termColumns(terms_) > 1)
            {
              model.e[ch] = x[coefIndex(ch, order_ + 1)];
            }
          }
        }

      private:
        size_t coefIndex(size_t channel, size_t entry) const
        {
          return order_ + channel * numerator_ + entry;
        }

        std::vector<ComplexT> buildPoles(const RealT* x) const
        {
          std::vector<ComplexT> poles(order_);
          size_t                q = 0;
          while (q < order_)
          {
            if (structure_[q].imag() == RealT{0})
            {
              poles[q]  = ComplexT{x[q], RealT{0}};
              q        += 1;
              continue;
            }
            poles[q]      = ComplexT{x[q], x[q + 1]};
            poles[q + 1]  = std::conj(poles[q]);
            q            += 2;
          }
          return poles;
        }

        ComplexT channelResidual(const RealT*                 x,
                                 size_t                       m,
                                 size_t                       ch,
                                 const std::vector<ComplexT>& phi) const
        {
          ComplexT value{RealT{0}, RealT{0}};
          for (size_t k = 0; k < order_; ++k)
          {
            value += x[coefIndex(ch, k)] * phi[k];
          }
          if (termColumns(terms_) > 0)
          {
            value += x[coefIndex(ch, order_)];
          }
          if (termColumns(terms_) > 1)
          {
            value += ComplexT{RealT{0}, samples_.omega[m]} * x[coefIndex(ch, order_ + 1)];
          }
          return value - samples_.response[m * channels_ + ch];
        }

        const SampledResponse<ScalarT, IdxT>& samples_;
        const std::vector<RealT>&             weights_;
        std::vector<ComplexT>                 structure_;
        RationalTerms                         terms_;
        RealT                                 margin_;
        size_t                                channels_;
        size_t                                order_;
        size_t                                numerator_;
        std::vector<RealT>                    start_;
      };

      /**
       * @brief Run the refinement program and keep the result only when
       *        it improves the weighted misfit.
       */
      template <typename scalar_type, typename index_type>
      int refineModel(
          const SampledResponse<scalar_type, index_type>&          samples,
          const std::vector<typename SampledResponse<
              scalar_type,
              index_type>::RealT>&                                 weights,
          RationalTerms                                            terms,
          typename SampledResponse<scalar_type, index_type>::RealT margin,
          RationalModel<scalar_type, index_type>&                  model)
      {
        using RealT   = typename SampledResponse<scalar_type,
                                                 index_type>::RealT;
        using SolverT = IpoptSolver<scalar_type, index_type>;

        RefinementProblem<scalar_type, index_type> problem(
            samples, weights, model, terms, margin);

        std::vector<RealT> start(
            static_cast<size_t>(problem.variableCount()), RealT{0});
        if (problem.startingPoint(start.data()) != 0)
        {
          return -2;
        }
        RealT before = RealT{0};
        if (problem.evaluateObjective(start.data(), before) != 0)
        {
          return -2;
        }

        typename SolverT::Parameters options;
        options.tolerance = REFINE_TOLERANCE;
        options.max_iterations =
            static_cast<index_type>(REFINE_MAX_ITERATIONS);
        options.hessian = HessianMode::LIMITED_MEMORY;

        SolverT            solver(&problem, options);
        std::vector<RealT> x = start;
        if (solver.solve(x) < 0)
        {
          return -2;
        }

        RealT after = RealT{0};
        if (problem.evaluateObjective(x.data(), after) != 0)
        {
          return -2;
        }
        if (after < before)
        {
          problem.unpack(x.data(), model);
        }
        return 0;
      }
    } // namespace Detail

    template <typename scalar_type, typename index_type>
    std::string VectorFitting<scalar_type, index_type>::Stats::report() const
    {
      std::ostringstream stream;
      stream << "VectorFitting: " << iterations.size() << " passes, "
             << (converged ? "converged" : "not converged")
             << ", rel rms " << final_rel_rms << ", order "
             << order_selected;
      return stream.str();
    }

    template <typename scalar_type, typename index_type>
    VectorFitting<scalar_type, index_type>::VectorFitting(
        const SampledResponse<ScalarT, IdxT>& samples)
      : samples_(samples)
    {
    }

    /**
     * @brief Fit a rational approximation to the sampled response.
     *
     * Orchestrates seeding, the relaxed relocation loop, coefficient
     * identification, deterministic restarts, the lowest-order search,
     * and the optional refinement stage. See README.md.
     */
    template <typename scalar_type, typename index_type>
    int VectorFitting<scalar_type, index_type>::fit(
        RationalModel<ScalarT, IdxT>& model, Parameters params)
    {
      const auto channels     = static_cast<size_t>(samples_.rows) * static_cast<size_t>(samples_.cols);
      const auto sample_count = samples_.omega.size();

      if (sample_count < 2 || samples_.response.size() != sample_count * channels)
      {
        return -1;
      }
      for (size_t m = 0; m < sample_count; ++m)
      {
        if (samples_.omega[m] <= RealT{0} || (m > 0 && samples_.omega[m] <= samples_.omega[m - 1]))
        {
          return -1;
        }
      }
      if ((params.weighting == Weighting::CUSTOM) != (params.weights.size() == sample_count))
      {
        return -1;
      }
      if ((params.seeding == PoleSeeding::USER) != !params.starting_poles.empty())
      {
        return -1;
      }
      if (params.seeding == PoleSeeding::USER && params.order_search.enabled)
      {
        return -1;
      }
      if (params.order_search.enabled
              ? (params.order_search.min_poles < 1
                 || params.order_search.max_poles < params.order_search.min_poles)
              : params.pole_count < 1)
      {
        return -1;
      }

      // Per-sample weights by the documented conventions.
      std::vector<RealT> sample_weights(sample_count, RealT{1});
      if (params.weighting == Weighting::CUSTOM)
      {
        sample_weights = params.weights;
      }
      else if (params.weighting != Weighting::UNIFORM)
      {
        RealT              peak = RealT{0};
        std::vector<RealT> norms(sample_count, RealT{0});
        for (size_t m = 0; m < sample_count; ++m)
        {
          RealT norm = RealT{0};
          for (size_t ch = 0; ch < channels; ++ch)
          {
            norm += std::norm(samples_.response[m * channels + ch]);
          }
          norms[m] = std::sqrt(norm);
          peak     = std::max(peak, norms[m]);
        }
        const RealT floor = RealT{WEIGHT_FLOOR_RATIO} * peak + RealT{SHIFT_NORM_FLOOR};
        for (size_t m = 0; m < sample_count; ++m)
        {
          const RealT magnitude = std::max(norms[m], floor);
          sample_weights[m] =
              params.weighting == Weighting::INVERSE_MAGNITUDE
                  ? RealT{1} / magnitude
                  : RealT{1} / (magnitude * magnitude);
        }
      }

      std::vector<IdxT> counts;
      if (params.order_search.enabled)
      {
        const IdxT step = params.seeding == PoleSeeding::LOG_REAL ? 1 : 2;
        for (IdxT count = params.order_search.min_poles;
             count <= params.order_search.max_poles;
             count += step)
        {
          counts.push_back(count);
        }
      }
      else
      {
        counts.push_back(params.pole_count);
      }

      Stats                        best_stats;
      RationalModel<ScalarT, IdxT> best_model;
      bool                         have_best = false;

      for (const auto count : counts)
      {
        std::vector<ComplexT> seed;
        if (params.seeding == PoleSeeding::USER)
        {
          seed = canonicalizePoles(params.starting_poles,
                                   params.stability_margin);
          if (seed.size() != static_cast<size_t>(count))
          {
            return -1;
          }
        }
        else
        {
          const RealT low  = samples_.omega.front();
          const RealT high = samples_.omega.back();
          const bool  complex_seed =
              params.seeding == PoleSeeding::LOG_COMPLEX;
          const auto pair_count =
              complex_seed ? static_cast<size_t>(count) / 2
                           : static_cast<size_t>(0);
          const auto real_count =
              static_cast<size_t>(count) - 2 * pair_count;

          std::vector<ComplexT> raw;
          const auto            points = pair_count + real_count;
          for (size_t k = 0; k < points; ++k)
          {
            const RealT fraction =
                points == 1 ? RealT{0}
                            : static_cast<RealT>(k) / static_cast<RealT>(points - 1);
            const RealT beta =
                std::exp(std::log(low) + fraction * (std::log(high) - std::log(low)));
            if (k < real_count)
            {
              raw.push_back(ComplexT{-beta, RealT{0}});
            }
            else
            {
              raw.push_back(
                  ComplexT{-RealT{SEED_DAMPING_RATIO} * beta, beta});
            }
          }
          seed = canonicalizePoles(raw, params.stability_margin);
        }

        Stats                        attempt_stats;
        RationalModel<ScalarT, IdxT> attempt_model;
        bool                         have_attempt = false;

        for (IdxT attempt = 0; attempt <= params.restarts.max_restarts;
             ++attempt)
        {
          if (attempt > 0)
          {
            // Deterministic bounded perturbation of the previous poles.
            const RealT phase = std::fmod(
                static_cast<RealT>(attempt) * RealT{RESTART_GOLDEN_FRACT},
                RealT{1});
            const RealT factor = RealT{0.7} + RealT{0.6} * phase;
            for (auto& pole : seed)
            {
              pole = ComplexT{pole.real() * factor, pole.imag() * factor};
            }
            seed = canonicalizePoles(seed, params.stability_margin);
          }

          Stats trial_stats;
          auto  poles = seed;
          for (IdxT pass = 0; pass < params.max_iterations; ++pass)
          {
            RealT     shift      = RealT{0};
            bool      degenerate = false;
            const int relocation =
                Detail::relocatePoles<ScalarT, IdxT>(samples_,
                                                     sample_weights,
                                                     params.terms,
                                                     params.stability_margin,
                                                     poles,
                                                     shift,
                                                     degenerate);
            if (relocation != 0)
            {
              return relocation;
            }
            trial_stats.iterations.push_back({pass, shift});
            if (degenerate)
            {
              // A null relocation is a fixed point without a convergence
              // certificate; iterating further cannot move the poles.
              break;
            }
            if (shift < params.pole_shift_tolerance)
            {
              trial_stats.converged = true;
              break;
            }
          }

          RationalModel<ScalarT, IdxT> trial_model;
          const int                    identification =
              Detail::identifyCoefficients<ScalarT, IdxT>(samples_,
                                                          sample_weights,
                                                          poles,
                                                          params.terms,
                                                          trial_model);
          if (identification != 0)
          {
            return identification;
          }

          Detail::computeStats<ScalarT, IdxT>(samples_, trial_model, trial_stats);
          trial_stats.restarts_used = attempt;

          if (!have_attempt || trial_stats.final_rel_rms < attempt_stats.final_rel_rms)
          {
            attempt_stats = trial_stats;
            attempt_model = trial_model;
            have_attempt  = true;
            seed          = poles;
          }
          if (attempt_stats.final_rel_rms <= params.restarts.trigger_rel_rms)
          {
            break;
          }
        }

        if (!have_best || attempt_stats.final_rel_rms < best_stats.final_rel_rms)
        {
          best_stats                = attempt_stats;
          best_model                = attempt_model;
          best_stats.order_selected = count;
          have_best                 = true;
        }
        if (params.order_search.enabled && best_stats.final_rel_rms <= params.order_search.target_rel_rms)
        {
          break;
        }
      }

      model  = best_model;
      stats_ = best_stats;

      if (params.refine)
      {
        const int refinement =
            Detail::refineModel<ScalarT, IdxT>(samples_,
                                               sample_weights,
                                               params.terms,
                                               params.stability_margin,
                                               model);
        if (refinement != 0)
        {
          return refinement;
        }
        Detail::computeStats<ScalarT, IdxT>(samples_, model, stats_);
      }

      if (model.validate(RealT{0}) != 0)
      {
        return -3;
      }
      if (params.order_search.enabled && stats_.final_rel_rms > params.order_search.target_rel_rms)
      {
        return 2;
      }
      return stats_.converged ? 0 : 1;
    }

    template class VectorFitting<double, long int>;
    template class VectorFitting<double, size_t>;
  } // namespace Optimization
} // namespace GridKit
