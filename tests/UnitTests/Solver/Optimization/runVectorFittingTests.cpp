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
  result += fit_refine_never_worsens_verdict();
  result += fit_rejects_invalid_order_ranges();
  return result.summary();
}
