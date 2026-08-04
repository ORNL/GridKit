/**
 * @file runModalDecompositionTests.cpp
 *
 * @brief Unit tests for the identity-preserving modal observation of
 *        the propagation matrix.
 *
 */

#include <cmath>
#include <complex>
#include <vector>

#include <GridKit/Model/EMT/Parameters/Modal/ModalDecomposition.hpp>
#include <GridKit/Testing/Testing.hpp>

#include <Eigen/Dense>

namespace
{
  using scalar_type = double;
  using index_type  = size_t;

  using ModalT =
      GridKit::EMT::Parameters::ModalDecomposition<scalar_type, index_type>;
  using ComplexT = std::complex<scalar_type>;

  bool close(scalar_type value, scalar_type reference, scalar_type tol)
  {
    return std::abs(value - reference) <= tol;
  }

  /// Row-major real and imaginary parts of a complex matrix, the form
  /// decompose() consumes.
  struct Parts
  {
    std::vector<scalar_type> alpha;
    std::vector<scalar_type> beta;
  };

  Parts partsOf(const Eigen::MatrixXcd& gamma)
  {
    const auto k = gamma.rows();
    Parts      parts;
    parts.alpha.resize(static_cast<size_t>(k * k));
    parts.beta.resize(static_cast<size_t>(k * k));
    for (Eigen::Index i = 0; i < k; ++i)
    {
      for (Eigen::Index j = 0; j < k; ++j)
      {
        const auto entry   = static_cast<size_t>(i * k + j);
        parts.alpha[entry] = gamma(i, j).real();
        parts.beta[entry]  = gamma(i, j).imag();
      }
    }
    return parts;
  }

  /// The published sample read back into Eigen form.
  struct Frame
  {
    Eigen::MatrixXcd tv;
    Eigen::MatrixXcd ti;
    Eigen::VectorXcd lambda;
    Eigen::VectorXcd h;
    Eigen::VectorXd  tau;
  };

  Frame frameOf(const ModalT& modes)
  {
    const auto  k      = static_cast<Eigen::Index>(modes.modeCount());
    const auto* values = modes.data();

    Frame frame{Eigen::MatrixXcd(k, k),
                Eigen::MatrixXcd(k, k),
                Eigen::VectorXcd(k),
                Eigen::VectorXcd(k),
                Eigen::VectorXd(k)};

    const auto tv_r = modes.tvReal().offset;
    const auto tv_i = modes.tvImag().offset;
    const auto ti_r = modes.tiReal().offset;
    const auto ti_i = modes.tiImag().offset;
    for (Eigen::Index i = 0; i < k; ++i)
    {
      for (Eigen::Index j = 0; j < k; ++j)
      {
        const auto entry = static_cast<size_t>(i * k + j);
        frame.tv(i, j)   = ComplexT{values[tv_r + entry],
                                  values[tv_i + entry]};
        frame.ti(i, j)   = ComplexT{values[ti_r + entry],
                                  values[ti_i + entry]};
      }
    }
    for (Eigen::Index m = 0; m < k; ++m)
    {
      const auto entry = static_cast<size_t>(m);
      frame.lambda(m)  = ComplexT{values[modes.alphaM().offset + entry],
                                 values[modes.betaM().offset + entry]};
      frame.h(m)       = ComplexT{values[modes.hReal().offset + entry],
                            values[modes.hImag().offset + entry]};
      frame.tau(m)     = values[modes.tau().offset + entry];
    }
    return frame;
  }

  scalar_type biorthogonalityDefect(const Frame& frame)
  {
    const auto k = frame.tv.rows();
    return (frame.ti.adjoint() * frame.tv
            - Eigen::MatrixXcd::Identity(k, k))
        .cwiseAbs()
        .maxCoeff();
  }

  /// A fixed, well-conditioned, deliberately nonorthogonal basis.
  Eigen::MatrixXcd testBasis()
  {
    Eigen::MatrixXcd basis(3, 3);
    basis << ComplexT{1.0, 0.1}, ComplexT{0.3, 0.0}, ComplexT{0.1, -0.2},
        ComplexT{0.2, -0.1}, ComplexT{1.0, 0.2}, ComplexT{0.3, 0.1},
        ComplexT{-0.1, 0.2}, ComplexT{0.2, -0.3}, ComplexT{1.0, 0.0};
    return basis;
  }

  /// Modes with distinct eigenvalues must reproduce the exact
  /// eigenpairs, machine-precision duals, and the closed-form modal
  /// responses.
  GridKit::Testing::TestOutcome decompose_matches_direct_eigensolve()
  {
    using GridKit::Testing::TestStatus;

    const scalar_type length = 1.0e5;
    const scalar_type omega  = 2.0 * M_PI * 50.0;

    Eigen::VectorXcd lambda(3);
    lambda << ComplexT{1.0e-5, 2.0e-4}, ComplexT{3.0e-5, 5.0e-4},
        ComplexT{6.0e-5, 9.0e-4};
    const Eigen::MatrixXcd basis = testBasis();
    const Eigen::MatrixXcd gamma =
        basis * lambda.asDiagonal() * basis.inverse();
    const auto parts = partsOf(gamma);

    ModalT modes(3, length);

    TestStatus success  = true;
    success            *= (modes.modeCount() == 3);
    success            *= (modes.decompose(omega, parts.alpha.data(), parts.beta.data())
                == 0);

    const auto frame  = frameOf(modes);
    success          *= (biorthogonalityDefect(frame) < 1.0e-12);

    for (Eigen::Index m = 0; m < 3; ++m)
    {
      // Eigen-relation residual certifies (lambda_m, tv_m) jointly.
      const scalar_type residual =
          (gamma * frame.tv.col(m) - frame.lambda(m) * frame.tv.col(m))
              .norm();
      success *= (residual < 1.0e-9);

      const ComplexT factor  = std::exp(-length * frame.lambda(m));
      success               *= close(frame.h(m).real(), factor.real(), 1.0e-10);
      success               *= close(frame.h(m).imag(), factor.imag(), 1.0e-10);
      success               *= close(frame.tau(m),
                       length * frame.lambda(m).imag() / omega,
                       1.0e-10);
    }

    // The canonical first frame ascends by the imaginary part.
    success *= (frame.lambda(0).imag() < frame.lambda(1).imag());
    success *= (frame.lambda(1).imag() < frame.lambda(2).imag());

    return success.report(__func__);
  }

  /// Identity must follow the eigenvectors through an exact eigenvalue
  /// crossing, including the sample where the values coincide and the
  /// matrix carries no eigenvector information at all.
  GridKit::Testing::TestOutcome decompose_carries_identity_through_crossing()
  {
    using GridKit::Testing::TestStatus;

    const scalar_type length = 1.0e5;
    const scalar_type omega  = 2.0 * M_PI * 50.0;

    Eigen::MatrixXcd basis(2, 2);
    basis << ComplexT{1.0, 0.0}, ComplexT{0.4, 0.1},
        ComplexT{0.3, -0.2}, ComplexT{1.0, 0.0};

    ModalT modes(2, length);

    TestStatus       success      = true;
    Eigen::Index     moving_label = 0;
    Eigen::MatrixXcd product_previous;
    for (int step = 0; step <= 20; ++step)
    {
      const scalar_type t = 0.05 * static_cast<scalar_type>(step);

      // The value curves cross at t = 0.5, which the grid hits exactly.
      const scalar_type moving = (1.0 + 2.0 * t) * 1.0e-5;
      const scalar_type still  = 2.0e-5;
      Eigen::VectorXcd  lambda(2);
      lambda << ComplexT{moving, 2.0e-4}, ComplexT{still, 2.0e-4};
      const Eigen::MatrixXcd gamma =
          basis * lambda.asDiagonal() * basis.inverse();
      const auto parts = partsOf(gamma);

      success          *= (modes.decompose(omega, parts.alpha.data(), parts.beta.data())
                  == 0);
      const auto frame  = frameOf(modes);

      success *= (biorthogonalityDefect(frame) < 1.0e-10);

      // The first sample decides which label owns the moving curve;
      // both traces must then follow their own curve straight through
      // the crossing, where the values coincide and either label is
      // trivially consistent.
      if (step == 0)
      {
        moving_label = std::abs(frame.lambda(0).real() - moving)
                               < std::abs(frame.lambda(1).real() - moving)
                           ? 0
                           : 1;
      }
      success *= close(frame.lambda(moving_label).real(), moving, 1.0e-9);
      success *= close(frame.lambda(1 - moving_label).real(), still, 1.0e-9);

      // The physical object stays exact and continuous even at the
      // degenerate sample, where the basis itself may legitimately
      // snap to the orthonormalized span.
      const scalar_type defect =
          (frame.tv * frame.lambda.asDiagonal() * frame.ti.adjoint()
           - gamma)
              .cwiseAbs()
              .maxCoeff();
      success *= (defect < 1.0e-10);

      const Eigen::MatrixXcd product =
          frame.tv * frame.h.asDiagonal() * frame.ti.adjoint();
      if (step > 0)
      {
        const scalar_type drift =
            (product - product_previous).cwiseAbs().maxCoeff();
        success *= (drift < 0.25);
      }
      product_previous = product;
    }

    return success.report(__func__);
  }

  /// A persistently degenerate pair must stay continuous while its
  /// invariant subspace rotates: the semisimple symmetric-tower case
  /// that eigenvector continuation cannot traverse at all.
  GridKit::Testing::TestOutcome decompose_is_continuous_at_exact_degeneracy()
  {
    using GridKit::Testing::TestStatus;

    const scalar_type length = 1.0e5;
    const scalar_type omega  = 2.0 * M_PI * 50.0;

    ModalT modes(3, length);

    TestStatus       success = true;
    Eigen::MatrixXcd tv_previous;
    for (int step = 0; step <= 20; ++step)
    {
      const scalar_type t = 0.05 * static_cast<scalar_type>(step);

      // Two equal eigenvalues whose shared eigenspace rotates with t.
      Eigen::MatrixXd rotation = Eigen::MatrixXd::Identity(3, 3);
      const auto      angle    = 0.4 * t;
      rotation(0, 0)           = std::cos(angle);
      rotation(0, 2)           = std::sin(angle);
      rotation(2, 0)           = -std::sin(angle);
      rotation(2, 2)           = std::cos(angle);

      Eigen::VectorXcd lambda(3);
      lambda << ComplexT{(2.0 + t) * 1.0e-5, 3.0e-4},
          ComplexT{(2.0 + t) * 1.0e-5, 3.0e-4},
          ComplexT{6.0e-5, 8.0e-4};
      const Eigen::MatrixXcd gamma = rotation.cast<ComplexT>()
                                     * lambda.asDiagonal()
                                     * rotation.transpose().cast<ComplexT>();
      const auto parts = partsOf(gamma);

      success          *= (modes.decompose(omega, parts.alpha.data(), parts.beta.data())
                  == 0);
      const auto frame  = frameOf(modes);

      success *= (biorthogonalityDefect(frame) < 1.0e-10);

      // Semisimple degeneracy: the reported frame must reproduce Gamma
      // exactly even though its columns are only a subspace choice.
      const scalar_type defect =
          (frame.tv * frame.lambda.asDiagonal() * frame.ti.adjoint()
           - gamma)
              .cwiseAbs()
              .maxCoeff();
      success *= (defect < 1.0e-12);

      if (step > 0)
      {
        const scalar_type drift =
            (frame.tv - tv_previous).cwiseAbs().maxCoeff();
        success *= (drift < 0.1);
      }
      tv_previous = frame.tv;
    }

    return success.report(__func__);
  }

  /// A sub-gap cluster commits a bounded reconstruction error: the
  /// propagation factor rebuilt from the cluster basis must match the
  /// exact matrix exponential within the documented bound.
  GridKit::Testing::TestOutcome decompose_bounds_cluster_reconstruction()
  {
    using GridKit::Testing::TestStatus;

    const scalar_type length = 1.0e5;
    const scalar_type omega  = 2.0 * M_PI * 50.0;

    Eigen::MatrixXcd basis(2, 2);
    basis << ComplexT{1.0, 0.0}, ComplexT{0.5, 0.2},
        ComplexT{-0.3, 0.1}, ComplexT{1.0, 0.0};

    const ComplexT   base{2.0e-5, 3.0e-4};
    Eigen::VectorXcd lambda(2);
    lambda << base, base * (1.0 + 3.0e-9);
    const Eigen::MatrixXcd gamma =
        basis * lambda.asDiagonal() * basis.inverse();
    const auto parts = partsOf(gamma);

    ModalT modes(2, length);

    TestStatus success = true;
    for (int sample = 0; sample < 2; ++sample)
    {
      // The second call exercises the cluster alignment path.
      success *= (modes.decompose(omega, parts.alpha.data(), parts.beta.data())
                  == 0);
    }
    const auto frame = frameOf(modes);

    success *= (biorthogonalityDefect(frame) < 1.0e-10);

    const Eigen::MatrixXcd exact =
        basis
        * (-length * lambda.array()).exp().matrix().asDiagonal()
        * basis.inverse();
    const Eigen::MatrixXcd rebuilt =
        frame.tv * frame.h.asDiagonal() * frame.ti.adjoint();

    const scalar_type defect =
        (rebuilt - exact).cwiseAbs().maxCoeff() / exact.cwiseAbs().maxCoeff();
    success *= (defect < 1.0e-6);

    return success.report(__func__);
  }

  /// The first frame is deterministic, reset() restores it, and the
  /// contract rejects invalid construction and samples.
  GridKit::Testing::TestOutcome decompose_first_sample_is_canonical()
  {
    using GridKit::Testing::TestStatus;

    const scalar_type length = 1.0e5;
    const scalar_type omega  = 2.0 * M_PI * 50.0;

    Eigen::VectorXcd lambda(3);
    lambda << ComplexT{3.0e-5, 7.0e-4}, ComplexT{1.0e-5, 2.0e-4},
        ComplexT{6.0e-5, 4.0e-4};
    const Eigen::MatrixXcd basis = testBasis();
    const Eigen::MatrixXcd gamma =
        basis * lambda.asDiagonal() * basis.inverse();
    const auto parts = partsOf(gamma);

    ModalT first(3, length);
    ModalT second(3, length);

    TestStatus success  = true;
    success            *= (first.decompose(omega, parts.alpha.data(), parts.beta.data())
                == 0);
    success            *= (second.decompose(omega, parts.alpha.data(), parts.beta.data())
                == 0);
    for (index_type entry = 0; entry < 4 * 9 + 5 * 3; ++entry)
    {
      success *= (first.data()[entry] == second.data()[entry]);
    }

    // Perturb the tracked frame, then reset: the canonical frame must
    // come back bit for bit.
    ModalT tracked(3, length);
    success *= (tracked.decompose(omega, parts.alpha.data(), parts.beta.data())
                == 0);
    success *= (tracked.decompose(2.0 * omega, parts.alpha.data(), parts.beta.data())
                == 0);
    tracked.reset();
    success *= (tracked.decompose(omega, parts.alpha.data(), parts.beta.data())
                == 0);
    for (index_type entry = 0; entry < 4 * 9 + 5 * 3; ++entry)
    {
      success *= (tracked.data()[entry] == first.data()[entry]);
    }

    // Invalid samples are rejected, not absorbed.
    success *= (first.decompose(0.0, parts.alpha.data(), parts.beta.data())
                == -1);

    bool constructor_throws = false;
    try
    {
      ModalT invalid(0, length);
    }
    catch (const std::exception&)
    {
      constructor_throws = true;
    }
    success *= constructor_throws;

    return success.report(__func__);
  }
} // namespace

int main()
{
  GridKit::Testing::TestingResults result;
  result += decompose_matches_direct_eigensolve();
  result += decompose_carries_identity_through_crossing();
  result += decompose_is_continuous_at_exact_degeneracy();
  result += decompose_bounds_cluster_reconstruction();
  result += decompose_first_sample_is_canonical();
  return result.summary();
}
