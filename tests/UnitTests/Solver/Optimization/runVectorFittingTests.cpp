/**
 * @file runVectorFittingTests.cpp
 *
 * @brief Unit tests for the relaxed vector fitting solver.
 *
 */

#include <cmath>
#include <complex>
#include <limits>
#include <vector>

#include <GridKit/Solver/Optimization/VectorFitting/VectorFitting.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace
{
  using scalar_type = double;
  using index_type  = size_t;

  using FitterT   = GridKit::Optimization::VectorFitting<scalar_type, index_type>;
  using ModelT    = GridKit::Optimization::RationalModel<scalar_type, index_type>;
  using ResponseT = GridKit::Optimization::SampledResponse<scalar_type, index_type>;
  using ComplexT  = ResponseT::ComplexT;

  /// The exact three-pole reference: one real pole and one conjugate
  /// pair with a constant term, evaluated analytically.
  struct Reference
  {
    ComplexT    real_pole{-50.0, 0.0};
    scalar_type real_residue{20.0};
    ComplexT    pair_pole{-30.0, 400.0};
    ComplexT    pair_residue{5.0, 2.0};
    scalar_type d{1.5};

    ComplexT evaluate(scalar_type omega) const
    {
      const ComplexT s{0.0, omega};
      return real_residue / (s - real_pole) + pair_residue / (s - pair_pole)
             + std::conj(pair_residue) / (s - std::conj(pair_pole)) + d;
    }
  };

  ResponseT makeSamples(const Reference& reference)
  {
    constexpr size_t      count = 60;
    constexpr scalar_type low   = 1.0;
    constexpr scalar_type high  = 1.0e4;

    ResponseT samples;
    samples.rows = 1;
    samples.cols = 1;
    samples.omega.resize(count);
    samples.response.resize(count);
    for (size_t m = 0; m < count; ++m)
    {
      const scalar_type fraction =
          static_cast<scalar_type>(m) / static_cast<scalar_type>(count - 1);
      samples.omega[m] =
          std::exp(std::log(low) + fraction * (std::log(high) - std::log(low)));
      samples.response[m] = reference.evaluate(samples.omega[m]);
    }
    return samples;
  }

  scalar_type nearestPoleDistance(const ModelT& model, ComplexT reference)
  {
    scalar_type nearest = std::numeric_limits<scalar_type>::infinity();
    for (const auto& pole : model.poles)
    {
      nearest = std::min(nearest,
                         std::abs(pole - reference) / std::abs(reference));
    }
    return nearest;
  }

  /// The reference extended to a full matrix response: shared poles,
  /// per-channel residues and constants varied deterministically.
  ResponseT makeMatrixSamples(const Reference& reference,
                              index_type       rows,
                              index_type       cols)
  {
    const auto base     = makeSamples(reference);
    const auto channels = static_cast<size_t>(rows) * static_cast<size_t>(cols);

    ResponseT samples;
    samples.rows  = rows;
    samples.cols  = cols;
    samples.omega = base.omega;
    samples.response.resize(base.omega.size() * channels);
    for (size_t m = 0; m < base.omega.size(); ++m)
    {
      for (size_t ch = 0; ch < channels; ++ch)
      {
        Reference element    = reference;
        element.real_residue = reference.real_residue
                               + 3.0 * static_cast<scalar_type>(ch);
        element.pair_residue = reference.pair_residue
                               + ComplexT{0.7, -0.4}
                                     * static_cast<scalar_type>(ch);
        element.d = reference.d + 0.25 * static_cast<scalar_type>(ch);

        samples.response[m * channels + ch] =
            element.evaluate(samples.omega[m]);
      }
    }
    return samples;
  }

  /// Exact rational data must converge with certificate: the relocation
  /// fixed point is reachable, the shift metric must see it, and the
  /// equilibrated programs must resolve the coefficients.
  GridKit::Testing::TestOutcome fit_recovers_exact_rational()
  {
    using GridKit::Testing::TestStatus;

    const Reference reference;
    const auto      samples = makeSamples(reference);

    FitterT fitter(samples);

    FitterT::Parameters params;
    params.pole_count = 3;
    params.terms      = GridKit::Optimization::RationalTerms::CONSTANT;

    ModelT    model;
    const int status = fitter.fit(model, params);

    TestStatus success  = true;
    success            *= (status == 0);
    success            *= (fitter.getStats().converged == true);
    success            *= (fitter.getStats().final_rel_rms < 1.0e-8);
    success            *= (model.poles.size() == 3);
    success            *= (nearestPoleDistance(model, reference.real_pole) < 1.0e-6);
    success            *= (nearestPoleDistance(model, reference.pair_pole) < 1.0e-6);
    success            *= (model.d.size() == 1);
    if (!model.d.empty())
    {
      success *= (std::abs(model.d[0] - reference.d) < 1.0e-6);
    }

    // The report surfaces the worst per-channel error so a collapsed
    // channel cannot hide inside the norm-weighted global figure.
    success *= (fitter.getStats().report().find("worst channel")
                != std::string::npos);

    return success.report(__func__);
  }

  /// A full matrix response with shared poles must be recovered
  /// exactly: the per-channel block elimination of the sigma program
  /// and the shared-basis identification both see nine channels here.
  GridKit::Testing::TestOutcome fit_recovers_exact_matrix_rational()
  {
    using GridKit::Testing::TestStatus;

    const Reference reference;
    const auto      samples = makeMatrixSamples(reference, 3, 3);

    FitterT fitter(samples);

    FitterT::Parameters params;
    params.pole_count = 3;
    params.terms      = GridKit::Optimization::RationalTerms::CONSTANT;

    ModelT    model;
    const int status = fitter.fit(model, params);

    TestStatus success  = true;
    success            *= (status == 0);
    success            *= (fitter.getStats().converged == true);
    success            *= (fitter.getStats().final_rel_rms < 1.0e-8);
    success            *= (model.poles.size() == 3);
    success            *= (nearestPoleDistance(model, reference.real_pole) < 1.0e-6);
    success            *= (nearestPoleDistance(model, reference.pair_pole) < 1.0e-6);
    success            *= (model.d.size() == 9);
    for (size_t ch = 0; ch < model.d.size(); ++ch)
    {
      const scalar_type expected =
          reference.d + 0.25 * static_cast<scalar_type>(ch);
      success *= (std::abs(model.d[ch] - expected) < 1.0e-6);
    }

    return success.report(__func__);
  }

  /// Refinement is accepted under the same unweighted metric that
  /// decides the verdict, so enabling it can never worsen the reported
  /// error or flip a met target into a failure.
  GridKit::Testing::TestOutcome fit_refine_never_worsens_verdict()
  {
    using GridKit::Testing::TestStatus;

    const Reference reference;
    const auto      samples = makeSamples(reference);

    FitterT::Parameters params;
    params.pole_count = 3;
    params.terms      = GridKit::Optimization::RationalTerms::CONSTANT;
    params.weighting  = GridKit::Optimization::Weighting::INVERSE_MAGNITUDE;

    FitterT   plain_fitter(samples);
    ModelT    plain_model;
    const int plain_status = plain_fitter.fit(plain_model, params);

    params.refine = true;
    FitterT   refined_fitter(samples);
    ModelT    refined_model;
    const int refined_status = refined_fitter.fit(refined_model, params);

    TestStatus success  = true;
    success            *= (plain_status >= 0);
    success            *= (refined_status >= 0);
    success            *= (refined_fitter.getStats().final_rel_rms
                <= plain_fitter.getStats().final_rel_rms * (1.0 + 1.0e-12));

    return success.report(__func__);
  }

  /// The reference matrix response with one channel replaced by
  /// deterministic noise at the numerical floor, the shape a monitored
  /// structural zero arrives in.
  ResponseT makeNoisyChannelSamples(const Reference& reference)
  {
    auto samples = makeMatrixSamples(reference, 2, 2);

    const auto channels = static_cast<size_t>(samples.rows)
                          * static_cast<size_t>(samples.cols);
    for (size_t m = 0; m < samples.omega.size(); ++m)
    {
      const auto phase = static_cast<scalar_type>(m);
      samples.response[m * channels + (channels - 1)] =
          ComplexT{1.0e-12 * std::cos(3.7 * phase),
                   1.0e-12 * std::sin(2.9 * phase)};
    }
    return samples;
  }

  /// Values below the relative floor are fit as exact zeros: the noise
  /// channel is flagged structural, its error measures leakage against
  /// the floor instead of a noise ratio, and the live channels still
  /// recover exactly. Without the floor the same channel poisons the
  /// per-channel metric.
  GridKit::Testing::TestOutcome fit_min_magnitude_declares_structural_zero()
  {
    using GridKit::Testing::TestStatus;

    const Reference reference;
    const auto      samples  = makeNoisyChannelSamples(reference);
    const auto      channels = static_cast<size_t>(samples.rows)
                          * static_cast<size_t>(samples.cols);

    FitterT::Parameters params;
    params.pole_count = 3;
    params.terms      = GridKit::Optimization::RationalTerms::CONSTANT;

    FitterT   raw_fitter(samples);
    ModelT    raw_model;
    const int raw_status = raw_fitter.fit(raw_model, params);

    params.min_mag = 1.0e-4;
    FitterT   cleaned_fitter(samples);
    ModelT    cleaned_model;
    const int cleaned_status = cleaned_fitter.fit(cleaned_model, params);

    const auto& raw     = raw_fitter.getStats();
    const auto& cleaned = cleaned_fitter.getStats();

    TestStatus success  = true;
    success            *= (raw_status >= 0);
    success            *= (cleaned_status == 0);

    // Without the floor the noise channel dominates the per-channel
    // metric and nothing is flagged.
    success *= (raw.channel_rel_rms[channels - 1] > 0.1);
    success *= (raw.report().find("structural zeros") == std::string::npos);

    // With the floor the channel is a declared zero: recovered exactly,
    // flagged, and reported.
    success *= (cleaned.final_rel_rms < 1.0e-8);
    success *= (cleaned.channel_structural_zero.size() == channels);
    for (size_t ch = 0; ch + 1 < channels; ++ch)
    {
      success *= (cleaned.channel_structural_zero[ch] == false);
    }
    success *= (cleaned.channel_structural_zero[channels - 1] == true);
    success *= (cleaned.channel_rel_rms[channels - 1] < 1.0e-8);
    success *= (cleaned.report().find("structural zeros 1")
                != std::string::npos);

    // The fitted model puts nothing where the target is zero.
    const auto rows = static_cast<index_type>(samples.rows);
    for (const auto omega : samples.omega)
    {
      success *= (std::abs(cleaned_model.evaluate(omega,
                                                  rows - 1,
                                                  rows - 1))
                  < 1.0e-10);
    }

    return success.report(__func__);
  }

  /// The relative floor must lie in [0, 1); anything else is a hard
  /// error, never a silent clamp.
  GridKit::Testing::TestOutcome fit_rejects_invalid_min_magnitude()
  {
    using GridKit::Testing::TestStatus;

    const Reference reference;
    const auto      samples = makeSamples(reference);

    FitterT fitter(samples);
    ModelT  model;

    FitterT::Parameters negative;
    negative.min_mag = -0.1;

    FitterT::Parameters saturated;
    saturated.min_mag = 1.0;

    FitterT::Parameters undefined;
    undefined.min_mag =
        std::numeric_limits<scalar_type>::quiet_NaN();

    TestStatus success  = true;
    success            *= (fitter.fit(model, negative) == -1);
    success            *= (fitter.fit(model, saturated) == -1);
    success            *= (fitter.fit(model, undefined) == -1);

    return success.report(__func__);
  }

  /// Invalid order ranges are hard errors, never a default-constructed
  /// model reported as success.
  GridKit::Testing::TestOutcome fit_rejects_invalid_order_ranges()
  {
    using GridKit::Testing::TestStatus;

    const Reference reference;
    const auto      samples = makeSamples(reference);

    FitterT fitter(samples);
    ModelT  model;

    FitterT::Parameters empty_range;
    empty_range.order_search.enabled   = true;
    empty_range.order_search.min_poles = 5;
    empty_range.order_search.max_poles = 2;

    FitterT::Parameters zero_min;
    zero_min.order_search.enabled   = true;
    zero_min.order_search.min_poles = 0;

    FitterT::Parameters zero_count;
    zero_count.pole_count = 0;

    TestStatus success  = true;
    success            *= (fitter.fit(model, empty_range) == -1);
    success            *= (fitter.fit(model, zero_min) == -1);
    success            *= (fitter.fit(model, zero_count) == -1);

    return success.report(__func__);
  }
} // namespace

int main()
{
  GridKit::Testing::TestingResults result;
  result += fit_recovers_exact_rational();
  result += fit_recovers_exact_matrix_rational();
  result += fit_refine_never_worsens_verdict();
  result += fit_min_magnitude_declares_structural_zero();
  result += fit_rejects_invalid_min_magnitude();
  result += fit_rejects_invalid_order_ranges();
  return result.summary();
}
