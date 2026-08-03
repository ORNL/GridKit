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
  constexpr double REFINE_TOLERANCE      = 1.0e-8;
  constexpr size_t REFINE_MAX_ITERATIONS = 300;
  constexpr double RESTART_GOLDEN_FRACT  = 0.6180339887;

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

  /// Sampled responses viewed as a samples-by-channels matrix.
  using ResponseView =
      Eigen::Map<const Eigen::Matrix<std::complex<double>,
                                     Eigen::Dynamic,
                                     Eigen::Dynamic,
                                     Eigen::RowMajor>>;

  /**
   * @brief Basis values for every sample: entry (m, k) holds basis
   *        function k evaluated at s = j omega_m.
   */
  Eigen::MatrixXcd
  basisTable(const Eigen::Ref<const Eigen::VectorXd>& omega,
             const std::vector<std::complex<double>>& poles)
  {
    Eigen::MatrixXcd                  table(omega.size(),
                           static_cast<Eigen::Index>(poles.size()));
    std::vector<std::complex<double>> phi(poles.size());
    for (Eigen::Index m = 0; m < omega.size(); ++m)
    {
      fillBasis(poles, omega(m), phi);
      table.row(m) =
          Eigen::Map<const Eigen::RowVectorXcd>(phi.data(), table.cols());
    }
    return table;
  }

  /**
   * @brief Weighted real-augmented regression matrix over the numerator
   *        columns: the top half carries sqrt(w_m) times the real part
   *        of each column at sample m, the bottom half the imaginary
   *        part. Pole columns come first, then the optional constant
   *        and linear-in-omega columns. Both estimation stages regress
   *        on this matrix.
   */
  Eigen::MatrixXd
  regressionMatrix(const Eigen::MatrixXcd&                  table,
                   const Eigen::Ref<const Eigen::VectorXd>& omega,
                   const Eigen::VectorXd&                   root_weight,
                   RationalTerms                            terms)
  {
    const auto samples = table.rows();
    const auto order   = table.cols();
    const auto columns =
        order + static_cast<Eigen::Index>(termColumns(terms));

    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(2 * samples, columns);
    matrix.topLeftCorner(samples, order) =
        root_weight.asDiagonal() * table.real();
    matrix.bottomLeftCorner(samples, order) =
        root_weight.asDiagonal() * table.imag();
    if (termColumns(terms) > 0)
    {
      matrix.col(order).head(samples) = root_weight;
    }
    if (termColumns(terms) > 1)
    {
      matrix.col(order + 1).tail(samples) = root_weight.cwiseProduct(omega);
    }
    return matrix;
  }

  /**
   * @brief Reciprocal column norms as equilibration scales; a column
   *        below the norm floor keeps unit scale.
   */
  Eigen::VectorXd columnScales(const Eigen::VectorXd& norms)
  {
    return norms.unaryExpr(
        [](double norm)
        { return norm > COLUMN_NORM_FLOOR ? 1.0 / norm : 1.0; });
  }

  /**
   * @brief Minimizer of ||matrix x|| subject to direction . x = target.
   *
   * The Householder rotation taking direction onto its leading
   * coordinate turns the constraint into a fixed leading component of
   * the rotated variables and leaves their tail as an unconstrained
   * least squares, so the equality is honored exactly without an
   * iterative solver.
   *
   * @return the minimizer, or an empty vector when the constraint
   *         direction degenerates
   */
  Eigen::VectorXd constrainedMinimizer(const Eigen::MatrixXd& matrix,
                                       const Eigen::VectorXd& direction,
                                       double                 target)
  {
    const Eigen::HouseholderQR<Eigen::MatrixXd> rotation(direction);

    const double leading = rotation.matrixQR()(0, 0);
    if (std::abs(leading) < COLUMN_NORM_FLOOR)
    {
      return {};
    }

    const Eigen::MatrixXd rotated  = matrix * rotation.householderQ();
    const auto            trailing = matrix.cols() - 1;

    Eigen::VectorXd z(matrix.cols());
    z(0)             = target / leading;
    z.tail(trailing) = rotated.rightCols(trailing)
                           .colPivHouseholderQr()
                           .solve(-z(0) * rotated.col(0));
    return rotation.householderQ() * z;
  }
} // namespace

namespace GridKit
{
  namespace Optimization
  {
    namespace Detail
    {
      /**
       * @brief Solve the relaxed sigma program by block elimination.
       *
       * The stage-A least squares is block arrow: the rows of a channel
       * touch only that channel's numerator coefficients and the shared
       * sigma coefficients, and the numerator basis is the same for
       * every channel. A thin QR per channel therefore eliminates the
       * numerator block exactly, the surviving triangles stack into a
       * reduced problem over sigma alone, and the relaxation equality
       * is absorbed by a Householder change of variables. The result is
       * the exact minimizer of the fast formulation of Deschrijver et
       * al.: cost grows linearly in the channel count, no normal
       * equations are formed, and no iterative solver runs. Only sigma
       * is returned, because relocation consumes nothing else.
       *
       * Columns are equilibrated first: numerator columns by their own
       * weighted norm, sigma columns by their norm over all channels,
       * since the stacked triangles share the sigma variables and
       * per-channel scales would make the blocks incommensurate.
       *
       * @return 0 on success with sigma of size order + 1; -2 when the
       *         samples cannot determine the requested order or the
       *         relaxation direction degenerates
       */
      template <typename scalar_type, typename index_type>
      int solveSigmaProgram(
          const SampledResponse<scalar_type, index_type>& samples,
          const std::vector<typename SampledResponse<
              scalar_type,
              index_type>::RealT>&                        sample_weights,
          const std::vector<typename SampledResponse<
              scalar_type,
              index_type>::ComplexT>&                     poles,
          RationalTerms                                   terms,
          std::vector<typename SampledResponse<
              scalar_type,
              index_type>::RealT>&                        sigma)
      {
        using RealT = typename SampledResponse<scalar_type,
                                               index_type>::RealT;

        const auto channels = static_cast<Eigen::Index>(samples.rows)
                              * static_cast<Eigen::Index>(samples.cols);
        const auto sample_count =
            static_cast<Eigen::Index>(samples.omega.size());
        const auto order = static_cast<Eigen::Index>(poles.size());
        const auto numerator =
            order + static_cast<Eigen::Index>(termColumns(terms));
        const auto sigma_count = order + 1;

        if (2 * sample_count < numerator + sigma_count)
        {
          return -2;
        }

        const Eigen::Map<const Eigen::VectorXd> omega(samples.omega.data(),
                                                      sample_count);
        const Eigen::Map<const Eigen::VectorXd> weights(
            sample_weights.data(), sample_count);
        const ResponseView    response(samples.response.data(), sample_count, channels);
        const Eigen::VectorXd root_weight = weights.cwiseSqrt();

        const Eigen::MatrixXcd table = basisTable(omega, poles);

        Eigen::MatrixXd basis =
            regressionMatrix(table, omega, root_weight, terms);
        basis = basis
                * columnScales(basis.colwise().norm().transpose())
                      .asDiagonal();

        // The weighted norm of sigma column k over every channel's rows
        // collapses to a sample sum against the per-sample response
        // energy, so the global scales come without touching a channel.
        const Eigen::VectorXd energy =
            weights.cwiseProduct(response.rowwise().squaredNorm());
        Eigen::VectorXd sigma_norms(sigma_count);
        sigma_norms.head(order) =
            (energy.transpose() * table.cwiseAbs2()).transpose().cwiseSqrt();
        sigma_norms(order)                 = std::sqrt(energy.sum());
        const Eigen::VectorXd sigma_scales = columnScales(sigma_norms);

        // The relaxation equality lives on the sigma block alone and
        // rides the same equilibration.
        Eigen::VectorXd relaxation(sigma_count);
        relaxation.head(order) = table.real().colwise().sum().transpose();
        relaxation(order)      = static_cast<RealT>(sample_count);
        relaxation             = relaxation.cwiseProduct(sigma_scales);

        // Eliminate each channel's numerator block. The QR must stay
        // unpivoted: column pivoting would mix numerator and sigma
        // columns and break the block split the elimination rests on.
        Eigen::MatrixXd stacked(2 * sample_count, numerator + sigma_count);
        stacked.leftCols(numerator) = basis;

        Eigen::MatrixXd  reduced(channels * sigma_count, sigma_count);
        Eigen::MatrixXcd modulated(sample_count, sigma_count);
        for (Eigen::Index ch = 0; ch < channels; ++ch)
        {
          modulated.leftCols(order) = response.col(ch).asDiagonal() * table;
          modulated.col(order)      = response.col(ch);

          auto block = stacked.rightCols(sigma_count);
          block.topRows(sample_count) =
              -(root_weight.asDiagonal() * modulated.real())
              * sigma_scales.asDiagonal();
          block.bottomRows(sample_count) =
              -(root_weight.asDiagonal() * modulated.imag())
              * sigma_scales.asDiagonal();

          const Eigen::HouseholderQR<Eigen::MatrixXd> qr(stacked);
          reduced.middleRows(ch * sigma_count, sigma_count) =
              qr.matrixQR()
                  .block(numerator, numerator, sigma_count, sigma_count)
                  .triangularView<Eigen::Upper>();
        }

        const Eigen::VectorXd minimizer = constrainedMinimizer(
            reduced, relaxation, static_cast<RealT>(sample_count));
        if (minimizer.size() == 0)
        {
          return -2;
        }

        const Eigen::VectorXd solution =
            minimizer.cwiseProduct(sigma_scales);
        sigma.assign(solution.data(), solution.data() + solution.size());
        return 0;
      }

      /**
       * @brief One relaxed relocation pass: solve the sigma program by
       *        block elimination, then relocate the poles to the zeros
       *        of sigma via the real companion eigenproblem.
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

        const auto order = poles.size();

        std::vector<RealT> sigma;
        const int          solved = solveSigmaProgram<scalar_type, index_type>(
            samples, sample_weights, poles, terms, sigma);
        if (solved != 0)
        {
          return solved;
        }

        RealT sigma_scale = sigma[order];
        if (std::abs(sigma_scale) < RealT{SIGMA_SCALE_FLOOR})
        {
          sigma_scale = sigma_scale < RealT{0} ? -RealT{SIGMA_SCALE_FLOOR}
                                               : RealT{SIGMA_SCALE_FLOOR};
        }

        RealT gamma_norm = RealT{0};
        RealT pole_scale = RealT{0};
        for (size_t k = 0; k < order; ++k)
        {
          gamma_norm += sigma[k] * sigma[k];
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
            gamma(row)           = sigma[q];
            q                   += 1;
            continue;
          }
          companion(row, row)          = poles[q].real();
          companion(row, row + 1)      = poles[q].imag();
          companion(row + 1, row)      = -poles[q].imag();
          companion(row + 1, row + 1)  = poles[q].real();
          b_vec(row)                   = 2.0;
          gamma(row)                   = sigma[q];
          gamma(row + 1)               = sigma[q + 1];
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
       *        frozen: every channel regresses on the same weighted
       *        basis, so it is factored once and back-substituted per
       *        channel.
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

        const auto                              rows = static_cast<Eigen::Index>(sample_count);
        const Eigen::Map<const Eigen::VectorXd> omega(samples.omega.data(),
                                                      rows);
        const ResponseView                      response(samples.response.data(), rows, static_cast<Eigen::Index>(channels));
        const Eigen::VectorXd                   root_weight =
            Eigen::Map<const Eigen::VectorXd>(sample_weights.data(), rows)
                .cwiseSqrt();

        const Eigen::MatrixXcd table = basisTable(omega, poles);

        Eigen::MatrixXd basis =
            regressionMatrix(table, omega, root_weight, terms);
        const Eigen::VectorXd scales =
            columnScales(basis.colwise().norm().transpose());
        basis = basis * scales.asDiagonal();

        const Eigen::ColPivHouseholderQR<Eigen::MatrixXd> factorization(
            basis);

        std::vector<RealT> solution;
        Eigen::VectorXd    rhs(2 * rows);
        for (size_t ch = 0; ch < channels; ++ch)
        {
          const auto column = response.col(static_cast<Eigen::Index>(ch));
          rhs.head(rows)    = root_weight.cwiseProduct(column.real());
          rhs.tail(rows)    = root_weight.cwiseProduct(column.imag());

          const Eigen::VectorXd coefficients =
              scales.cwiseProduct(factorization.solve(rhs));
          solution.assign(coefficients.data(),
                          coefficients.data() + coefficients.size());

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
       * @brief Run the refinement program and unpack its candidate.
       *
       * Acceptance is not decided here: the caller compares the
       * candidate against the pre-refinement model under the same
       * unweighted metric that decides the fit verdict, so refinement
       * can never flip a met target into a failure.
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

        typename SolverT::Parameters options;
        options.tolerance = REFINE_TOLERANCE;
        options.max_iterations =
            static_cast<index_type>(REFINE_MAX_ITERATIONS);
        options.hessian = HessianMode::LIMITED_MEMORY;

        SolverT            solver(&problem, options);
        std::vector<RealT> x;
        if (solver.solve(x) < 0)
        {
          return -2;
        }

        problem.unpack(x.data(), model);
        return 0;
      }
    } // namespace Detail

    template <typename scalar_type, typename index_type>
    std::string VectorFitting<scalar_type, index_type>::Stats::report() const
    {
      // The global rel rms is norm-weighted, so a channel whose response
      // has collapsed can be arbitrarily wrong without moving it; the
      // worst per-channel figure is the one that sees such a channel.
      RealT worst = RealT{0};
      for (const auto value : channel_rel_rms)
      {
        worst = std::max(worst, value);
      }

      std::ostringstream stream;
      stream << "VectorFitting: " << iterations.size() << " passes, "
             << (converged ? "converged" : "not converged")
             << ", rel rms " << final_rel_rms << ", worst channel " << worst
             << ", order " << order_selected;
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
      if (params.order_search.enabled
          && (params.order_search.plateau_improvement < RealT{0}
              || params.order_search.plateau_improvement >= RealT{1}
              || (params.order_search.plateau_improvement > RealT{0}
                  && params.order_search.plateau_passes < 1)))
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

      RealT plateau_reference = std::numeric_limits<RealT>::infinity();
      IdxT  plateau_stalls    = 0;

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
        if (params.order_search.enabled
            && params.order_search.plateau_improvement > RealT{0})
        {
          if (attempt_stats.final_rel_rms
              < (RealT{1} - params.order_search.plateau_improvement)
                    * plateau_reference)
          {
            plateau_reference = attempt_stats.final_rel_rms;
            plateau_stalls    = 0;
          }
          else if (++plateau_stalls >= params.order_search.plateau_passes)
          {
            // Additional poles have stopped paying for themselves; the
            // best model and the order-search verdict stand as they
            // are.
            break;
          }
        }
      }

      model  = best_model;
      stats_ = best_stats;

      if (params.refine)
      {
        RationalModel<ScalarT, IdxT> refined = model;
        const int                    refinement =
            Detail::refineModel<ScalarT, IdxT>(samples_,
                                               sample_weights,
                                               params.terms,
                                               params.stability_margin,
                                               refined);
        if (refinement != 0)
        {
          return refinement;
        }

        // Accept the candidate only under the verdict metric itself.
        Stats refined_stats = stats_;
        Detail::computeStats<ScalarT, IdxT>(samples_, refined, refined_stats);
        if (refined_stats.final_rel_rms < stats_.final_rel_rms)
        {
          model  = refined;
          stats_ = refined_stats;
        }
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
