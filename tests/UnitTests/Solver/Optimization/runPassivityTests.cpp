/**
 * @file runPassivityTests.cpp
 *
 * @brief Unit tests for the rational-model passivity assessment.
 *
 */

#include <cmath>
#include <complex>
#include <limits>

#include <GridKit/Solver/Optimization/Rational/Passivity.hpp>
#include <GridKit/Solver/Optimization/Rational/RationalModel.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace
{
  using scalar_type = double;
  using index_type  = size_t;

  using ModelT   = GridKit::Optimization::RationalModel<scalar_type, index_type>;
  using ReportT  = GridKit::Optimization::PassivityReport<scalar_type, index_type>;
  using ComplexT = std::complex<scalar_type>;

  constexpr scalar_type band_low  = 1.0;
  constexpr scalar_type band_high = 1.0e3;

  ModelT makeScalarModel(ComplexT    pole,
                         ComplexT    residue,
                         scalar_type constant)
  {
    ModelT model;
    model.rows  = 1;
    model.cols  = 1;
    model.poles = {pole};
    model.residues.assign(1, residue);
    model.d.assign(1, constant);
    return model;
  }

  /// A positive-real RC admittance is certified passive and stable.
  GridKit::Testing::TestOutcome passivity_accepts_passive_model()
  {
    using GridKit::Testing::TestStatus;

    const auto model =
        makeScalarModel(ComplexT{-1.0, 0.0}, ComplexT{1.0, 0.0}, 0.5);

    ReportT   report;
    const int status = GridKit::Optimization::assessPassivity(model,
                                                              report,
                                                              band_low,
                                                              band_high);

    TestStatus success  = true;
    success            *= (status == 0);
    success            *= (report.stable == true);
    success            *= (report.passive == true);
    success            *= report.violations.empty();

    return success.report(__func__);
  }

  /// Y(0) = d + r/(-p) = 0.4 - 0.5 < 0 while Re Y(j omega) crosses zero
  /// at omega = 0.5: the classic DC violation the screen must anchor at
  /// zero exactly, with its upper edge located by bisection.
  GridKit::Testing::TestOutcome passivity_locates_dc_violation()
  {
    using GridKit::Testing::TestStatus;

    const auto model =
        makeScalarModel(ComplexT{-1.0, 0.0}, ComplexT{-0.5, 0.0}, 0.4);

    ReportT   report;
    const int status = GridKit::Optimization::assessPassivity(model,
                                                              report,
                                                              band_low,
                                                              band_high);

    TestStatus success  = true;
    success            *= (status == 0);
    success            *= (report.stable == true);
    success            *= (report.passive == false);
    success            *= (report.violations.size() == 1);
    if (report.violations.size() == 1)
    {
      success *= (report.violations[0].omega_start == 0.0);
      success *= (std::abs(report.violations[0].omega_end - 0.5) < 1.0e-6);
    }

    return success.report(__func__);
  }

  /// An indefinite constant term violates passivity at every frequency
  /// beyond the screen; the report extends that band to infinity.
  GridKit::Testing::TestOutcome passivity_flags_indefinite_constant_term()
  {
    using GridKit::Testing::TestStatus;

    ModelT model;
    model.rows  = 2;
    model.cols  = 2;
    model.poles = {ComplexT{-10.0, 0.0}};
    model.residues.assign(4, ComplexT{0.0, 0.0});
    model.d = {1.0, 0.0, 0.0, -0.1};

    ReportT   report;
    const int status = GridKit::Optimization::assessPassivity(model,
                                                              report,
                                                              band_low,
                                                              band_high);

    TestStatus success  = true;
    success            *= (status == 0);
    success            *= (report.passive == false);
    success            *= !report.violations.empty();
    if (!report.violations.empty())
    {
      success *= std::isinf(report.violations.back().omega_end);
    }

    return success.report(__func__);
  }

  /// An unstable model is never passive; the verdict gates on stability
  /// instead of screening a meaningless conductance floor.
  GridKit::Testing::TestOutcome passivity_gates_on_stability()
  {
    using GridKit::Testing::TestStatus;

    const auto model =
        makeScalarModel(ComplexT{1.0, 0.0}, ComplexT{1.0, 0.0}, 0.5);

    ReportT   report;
    const int status = GridKit::Optimization::assessPassivity(model,
                                                              report,
                                                              band_low,
                                                              band_high);

    TestStatus success  = true;
    success            *= (status == 0);
    success            *= (report.stable == false);
    success            *= (report.passive == false);

    return success.report(__func__);
  }

  /// Degenerate bands are invalid inputs, not silent verdicts.
  GridKit::Testing::TestOutcome passivity_rejects_invalid_band()
  {
    using GridKit::Testing::TestStatus;

    const auto model =
        makeScalarModel(ComplexT{-1.0, 0.0}, ComplexT{1.0, 0.0}, 0.5);

    ReportT    report;
    TestStatus success  = true;
    success            *= (GridKit::Optimization::assessPassivity(model, report, 0.0, band_high) == -1);
    success            *= (GridKit::Optimization::assessPassivity(model, report, band_high, band_low) == -1);

    return success.report(__func__);
  }
} // namespace

int main()
{
  GridKit::Testing::TestingResults result;
  result += passivity_accepts_passive_model();
  result += passivity_locates_dc_violation();
  result += passivity_flags_indefinite_constant_term();
  result += passivity_gates_on_stability();
  result += passivity_rejects_invalid_band();
  return result.summary();
}
