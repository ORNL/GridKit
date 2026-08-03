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

  /// A response with a deterministic noise floor keeps every target
  /// unreachable: the plateau stop must end the ladder early with the
  /// same verdict an exhausted search would carry.
  GridKit::Testing::TestOutcome fit_order_search_plateau_stops_early()
  {
    using GridKit::Testing::TestStatus;

    const Reference reference;
    auto            samples = makeSamples(reference);
    for (size_t m = 0; m < samples.response.size(); ++m)
    {
      const scalar_type ripple =
          1.0e-3 * std::sin(29.0 * static_cast<scalar_type>(m) + 1.0);
      samples.response[m] *= (1.0 + ripple);
    }

    FitterT::Parameters params;
    params.terms                       = GridKit::Optimization::RationalTerms::CONSTANT;
    params.order_search.enabled        = true;
    params.order_search.min_poles      = 4;
    params.order_search.max_poles      = 12;
    params.order_search.target_rel_rms = 1.0e-10;

    FitterT   exhaustive_fitter(samples);
    ModelT    exhaustive_model;
    const int exhaustive_status =
        exhaustive_fitter.fit(exhaustive_model, params);

    params.order_search.plateau_improvement = 0.3;
    params.order_search.plateau_passes      = 2;

    FitterT   plateau_fitter(samples);
    ModelT    plateau_model;
    const int plateau_status = plateau_fitter.fit(plateau_model, params);

    TestStatus success  = true;
    success            *= (exhaustive_status == 2);
    success            *= (plateau_status == 2);
    success            *= (plateau_fitter.getStats().order_selected
                < params.order_search.max_poles);
    success            *= (plateau_fitter.getStats().final_rel_rms < 0.1);
    success            *= (plateau_fitter.getStats().final_rel_rms > 1.0e-6);

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

    FitterT::Parameters excessive_plateau;
    excessive_plateau.order_search.enabled             = true;
    excessive_plateau.order_search.plateau_improvement = 1.5;

    FitterT::Parameters zero_plateau_passes;
    zero_plateau_passes.order_search.enabled             = true;
    zero_plateau_passes.order_search.plateau_improvement = 0.1;
    zero_plateau_passes.order_search.plateau_passes      = 0;

    TestStatus success  = true;
    success            *= (fitter.fit(model, empty_range) == -1);
    success            *= (fitter.fit(model, zero_min) == -1);
    success            *= (fitter.fit(model, zero_count) == -1);
    success            *= (fitter.fit(model, excessive_plateau) == -1);
    success            *= (fitter.fit(model, zero_plateau_passes) == -1);

    return success.report(__func__);
  }
} // namespace

int main()
{
  GridKit::Testing::TestingResults result;
  result += fit_recovers_exact_rational();
  result += fit_recovers_exact_matrix_rational();
  result += fit_order_search_plateau_stops_early();
  result += fit_refine_never_worsens_verdict();
  result += fit_rejects_invalid_order_ranges();
  return result.summary();
}
