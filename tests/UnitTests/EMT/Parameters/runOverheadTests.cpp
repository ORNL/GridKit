#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include <GridKit/Model/EMT/Constants.hpp>
#include <GridKit/Model/EMT/Parameters/Overhead.hpp>
#include <GridKit/Testing/Testing.hpp>

#include <Eigen/Dense>
#include <Eigen/LU>

namespace
{
  using scalar_type = double;
  using index_type  = size_t;

  constexpr scalar_type pi   = GridKit::EMT::Constants::pi<scalar_type>();
  constexpr scalar_type mu0  = GridKit::EMT::Constants::mu0<scalar_type>();
  constexpr scalar_type eps0 = GridKit::EMT::Constants::epsilon0<scalar_type>();

  bool close(scalar_type value, scalar_type ref, scalar_type tol = 1.0e-10)
  {
    return std::abs(value - ref) <= tol * (1.0 + std::abs(ref));
  }

  using Complex2x2 = std::array<std::complex<scalar_type>, 4>;

  Complex2x2 multiply2(const Complex2x2& A, const Complex2x2& B)
  {
    Complex2x2 product;
    for (index_type i = 0; i < 2; ++i)
    {
      for (index_type j = 0; j < 2; ++j)
      {
        product[i * 2 + j] = A[i * 2] * B[j] + A[i * 2 + 1] * B[2 + j];
      }
    }
    return product;
  }

  scalar_type maxAbs2(const Complex2x2& A)
  {
    scalar_type value = 0.0;
    for (const auto& entry : A)
    {
      value = std::max(value, std::abs(entry));
    }
    return value;
  }

  GridKit::EMT::Parameters::OverheadData<scalar_type, index_type> makeData()
  {
    GridKit::EMT::Parameters::OverheadData<scalar_type, index_type> data;
    data.path.length            = 100000.0;
    data.tower.K                = 2;
    data.tower.span             = 300.0;
    data.tower.height           = {10.0, 10.0};
    data.tower.position         = {0.0, 30.0};
    data.conductor.K            = 2;
    data.conductor.radius       = {0.01, 0.02};
    data.conductor.inner_radius = {0.0, 0.0};
    data.conductor.sigma        = {3.7e7, 5.8e7};
    data.conductor.mu           = {mu0, mu0};
    data.conductor.weight       = {11.0, 12.0};
    data.conductor.phase        = {"a", "b"};
    data.carson.earth_sigma     = 1.0e-6;
    data.carson.earth_eps       = 2.5 * eps0;
    return data;
  }

  GridKit::EMT::Parameters::OverheadData<scalar_type, index_type> makeScalarData()
  {
    auto data                   = makeData();
    data.tower.K                = 1;
    data.tower.height           = {10.0};
    data.tower.position         = {0.0};
    data.conductor.K            = 1;
    data.conductor.radius       = {0.01};
    data.conductor.inner_radius = {0.0};
    data.conductor.sigma        = {3.7e7};
    data.conductor.mu           = {mu0};
    data.conductor.weight       = {11.0};
    data.conductor.phase        = {"a"};
    return data;
  }

  /// Two twin bundles and one grounded shield: five conductors that
  /// reduce to two electrical terminals. The bundles are asymmetric so
  /// the reduced matrices have no accidental structure.
  GridKit::EMT::Parameters::OverheadData<scalar_type, index_type> makeBundledData()
  {
    GridKit::EMT::Parameters::OverheadData<scalar_type, index_type> data;
    data.path.length            = 100000.0;
    data.tower.K                = 5;
    data.tower.span             = 300.0;
    data.tower.height           = {20.0, 20.0, 22.0, 22.0, 30.0};
    data.tower.position         = {-10.0, -9.6, 10.0, 10.4, 0.0};
    data.conductor.K            = 5;
    data.conductor.radius       = {0.012, 0.012, 0.014, 0.014, 0.006};
    data.conductor.inner_radius = {0.0, 0.0, 0.0, 0.0, 0.0};
    data.conductor.sigma        = {3.7e7, 3.7e7, 3.7e7, 3.7e7, 5.0e6};
    data.conductor.mu           = {mu0, mu0, mu0, mu0, mu0};
    data.conductor.weight       = {11.0, 11.0, 12.0, 12.0, 4.0};
    data.conductor.phase        = {"a", "a", "b", "b", "g"};
    data.carson.earth_sigma     = 1.0e-2;
    data.carson.earth_eps       = 10.0 * eps0;
    return data;
  }

  GridKit::EMT::Parameters::OverheadData<scalar_type, index_type> makeThreeConductorData()
  {
    auto data                   = makeData();
    data.tower.K                = 3;
    data.tower.height           = {21.5, 25.0, 23.0};
    data.tower.position         = {-8.0, 0.4, 8.8};
    data.conductor.K            = 3;
    data.conductor.radius       = {0.0143, 0.0143, 0.0143};
    data.conductor.inner_radius = {0.0, 0.0, 0.0};
    data.conductor.sigma        = {3.5e7, 3.5e7, 3.5e7};
    data.conductor.mu           = {mu0, mu0, mu0};
    data.conductor.weight       = {12.5, 12.5, 12.5};
    data.conductor.phase        = {"a", "b", "c"};
    data.carson.earth_sigma     = 0.01;
    data.carson.earth_eps       = 10.0 * eps0;
    return data;
  }

  scalar_type determinant3(const std::array<scalar_type, 9>& A)
  {
    return A[0] * (A[4] * A[8] - A[5] * A[7])
           - A[1] * (A[3] * A[8] - A[5] * A[6])
           + A[2] * (A[3] * A[7] - A[4] * A[6]);
  }

  scalar_type csrValue(GridKit::LinearAlgebra::CsrMatrix<scalar_type, index_type>* csr,
                       index_type                                                  row,
                       index_type                                                  col)
  {
    index_type*  rowptr = csr->getRowData();
    index_type*  colind = csr->getColData();
    scalar_type* vals   = csr->getValues();
    for (index_type p = rowptr[row]; p < rowptr[row + 1]; ++p)
    {
      if (colind[p] == col)
      {
        return vals[p];
      }
    }
    return 0.0;
  }

  /// Whether the sparsity pattern stores an entry at (row, col); a
  /// stored zero and an absent entry are different structural claims.
  bool csrHasEntry(GridKit::LinearAlgebra::CsrMatrix<scalar_type, index_type>* csr,
                   index_type                                                  row,
                   index_type                                                  col)
  {
    index_type* rowptr = csr->getRowData();
    index_type* colind = csr->getColData();
    for (index_type p = rowptr[row]; p < rowptr[row + 1]; ++p)
    {
      if (colind[p] == col)
      {
        return true;
      }
    }
    return false;
  }

  using Overhead =
      GridKit::EMT::Parameters::Overhead<scalar_type, index_type>;

  scalar_type finiteDifferenceJacobianValue(Overhead&   model,
                                            index_type  row,
                                            index_type  col,
                                            scalar_type alpha)
  {
    auto* y  = model.y().getData();
    auto* yp = model.yp().getData();
    auto* f  = model.getResidual().getData();

    const scalar_type y0  = y[col];
    const scalar_type yp0 = yp[col];
    const scalar_type h   = std::sqrt(std::numeric_limits<scalar_type>::epsilon())
                          * (1.0 + std::max(std::abs(y0), std::abs(yp0)));

    auto residualAt = [&](scalar_type direction)
    {
      y[col]  = y0 + direction * h;
      yp[col] = yp0 + direction * alpha * h;
      model.evaluateResidual();
      return f[row];
    };

    const scalar_type fp = residualAt(1.0);
    const scalar_type fm = residualAt(-1.0);

    y[col]  = y0;
    yp[col] = yp0;
    model.evaluateResidual();

    return (fp - fm) / (2.0 * h);
  }

  /// The checked coupling must exist in the sparsity pattern: csrValue
  /// reads an absent entry as zero, and several reference derivatives
  /// here sit below the value tolerance, so without the structural
  /// requirement a Jacobian missing the entry entirely would still pass.
  bool jacobianMatchesFiniteDifference(Overhead&   model,
                                       index_type  row,
                                       index_type  col,
                                       scalar_type alpha,
                                       scalar_type tol = 1.0e-5)
  {
    model.evaluateJacobian();
    if (!csrHasEntry(model.getCsrJacobian(), row, col))
    {
      return false;
    }
    const scalar_type actual = csrValue(model.getCsrJacobian(), row, col);
    const scalar_type ref    = finiteDifferenceJacobianValue(model, row, col, alpha);
    return std::isfinite(actual) && close(actual, ref, tol);
  }

  struct CarsonReference
  {
    scalar_type R{0.0};
    scalar_type L{0.0};
  };

  CarsonReference carsonReference(const GridKit::EMT::Parameters::Overhead<scalar_type, index_type>& model,
                                  index_type                                                         i,
                                  index_type                                                         j)
  {
    const auto        gammaR = model.carson().gammaR();
    const auto        gammaI = model.carson().gammaI();
    const scalar_type gr     = model.y().getData()[gammaR.offset];
    const scalar_type gi     = model.y().getData()[gammaI.offset];
    const scalar_type d      = model.tower().horizontalDistance(i, j);
    const scalar_type s      = model.tower().height(i) + model.tower().height(j);
    const scalar_type dp     = model.tower().imageDistance(i, j);
    const scalar_type gg     = gr * gr + gi * gi;
    const scalar_type p_r    = gr / gg;
    const scalar_type p_i    = -gi / gg;
    const scalar_type h_r    = s + 2.0 * p_r;
    const scalar_type h_i    = 2.0 * p_i;
    const scalar_type A      = d * d + h_r * h_r - h_i * h_i;
    const scalar_type B      = 2.0 * h_r * h_i;
    const scalar_type M      = A * A + B * B;
    const scalar_type theta  = std::atan2(B, A);
    const scalar_type lambda = 0.25 * std::log(M) - std::log(dp);

    return {-model.omega() * mu0 * theta / (4.0 * pi), mu0 * lambda / (2.0 * pi)};
  }

  GridKit::Testing::TestOutcome tower_model()
  {
    using GridKit::Testing::TestStatus;

    auto                                                         data = makeData();
    GridKit::EMT::Parameters::Conductor<scalar_type, index_type> conductor(data.conductor);
    conductor.initialize();
    GridKit::EMT::Parameters::Tower<scalar_type, index_type> tower(data.tower, conductor);
    tower.initialize();

    const scalar_type image_distance = std::sqrt(30.0 * 30.0 + 20.0 * 20.0);

    TestStatus success  = true;
    success            *= (tower.size() == 0);
    success            *= (tower.conductorCount() == 2);
    success            *= close(tower.height(0), 10.0);
    success            *= close(tower.position(1), 30.0);
    success            *= close(tower.distance(0, 1), 30.0);
    success            *= close(tower.distance(1, 0), 30.0);
    success            *= close(tower.imageDistance(0, 1), image_distance);
    success            *= close(tower.imageDistance(1, 0), image_distance);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome geometric_inductance_model()
  {
    using GridKit::Testing::TestStatus;

    GridKit::EMT::Parameters::Overhead<scalar_type, index_type> model(makeData());
    model.updateTime(2.0 * pi * 50.0, 1.0);
    model.allocate();

    const auto        Lgeo           = model.geometricInductance().Lgeo();
    const scalar_type image_distance = std::sqrt(30.0 * 30.0 + 20.0 * 20.0);
    const scalar_type L00_ref        = mu0 * std::log(2.0 * 10.0 / 0.01) / (2.0 * pi);
    const scalar_type L01_ref        = mu0 * std::log(image_distance / 30.0) / (2.0 * pi);
    const scalar_type L11_ref        = mu0 * std::log(2.0 * 10.0 / 0.02) / (2.0 * pi);

    TestStatus success  = true;
    success            *= (Lgeo.rows == 2);
    success            *= (Lgeo.cols == 2);
    success            *= (Lgeo.size() == 4);
    success            *= close(model.y().getData()[Lgeo.offset], L00_ref);
    success            *= close(model.y().getData()[Lgeo.offset + 1], L01_ref);
    success            *= close(model.y().getData()[Lgeo.offset + 2], L01_ref);
    success            *= close(model.y().getData()[Lgeo.offset + 3], L11_ref);

    model.evaluateResidual();
    for (index_type i = 0; i < Lgeo.size(); ++i)
    {
      success *= close(model.getResidual().getData()[Lgeo.offset + i], 0.0);
    }

    success *= jacobianMatchesFiniteDifference(model, Lgeo.offset, Lgeo.offset, 1.0, 1.0e-7);

    auto sagged_data          = makeData();
    sagged_data.tower.tension = {200000.0, 120000.0};

    GridKit::EMT::Parameters::Overhead<scalar_type, index_type> sagged_model(sagged_data);
    sagged_model.updateTime(2.0 * pi * 50.0, 1.0);
    sagged_model.allocate();

    const auto        Lgeo_sagged = sagged_model.geometricInductance().Lgeo();
    const scalar_type L00_sagged_ref =
        mu0 * std::log(2.0 * sagged_model.tower().height(0) / sagged_data.conductor.radius[0])
        / (2.0 * pi);
    const scalar_type L01_sagged_ref =
        mu0 * std::log(sagged_model.tower().imageDistance(0, 1) / sagged_model.tower().distance(0, 1))
        / (2.0 * pi);
    success *= close(sagged_model.y().getData()[Lgeo_sagged.offset], L00_sagged_ref);
    success *= close(sagged_model.y().getData()[Lgeo_sagged.offset + 1], L01_sagged_ref);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome carson_model()
  {
    using GridKit::Testing::TestStatus;

    const scalar_type omega = 2.0 * pi * 50.0;
    auto              data  = makeData();

    GridKit::EMT::Parameters::Overhead<scalar_type, index_type> model(data);
    model.updateTime(omega, 1.0);
    model.allocate();

    const auto        eta         = model.carson().eta();
    const auto        gammaR      = model.carson().gammaR();
    const auto        gammaI      = model.carson().gammaI();
    const auto        Rcarson     = model.carson().Rcarson();
    const auto        Lcarson     = model.carson().Lcarson();
    const scalar_type a           = omega * omega * mu0 * data.carson.earth_eps;
    const scalar_type b           = omega * mu0 * data.carson.earth_sigma;
    const scalar_type eta_ref     = std::sqrt(a * a + b * b);
    const scalar_type gamma_r_ref = std::sqrt(0.5 * (eta_ref - a));
    const scalar_type gamma_i_ref = std::sqrt(0.5 * (eta_ref + a));

    TestStatus success  = true;
    success            *= (model.carson().size() == 11);
    success            *= close(model.y().getData()[eta.offset], eta_ref);
    success            *= close(model.y().getData()[gammaR.offset], gamma_r_ref);
    success            *= close(model.y().getData()[gammaI.offset], gamma_i_ref);
    success            *= (Rcarson.rows == 2);
    success            *= (Rcarson.cols == 2);
    success            *= (Lcarson.rows == 2);
    success            *= (Lcarson.cols == 2);

    for (index_type i = 0; i < Rcarson.rows; ++i)
    {
      for (index_type j = 0; j < Rcarson.cols; ++j)
      {
        const CarsonReference ref  = carsonReference(model, i, j);
        success                   *= close(model.y().getData()[Rcarson.offset + i * Rcarson.cols + j], ref.R, 1.0e-12);
        success                   *= close(model.y().getData()[Lcarson.offset + i * Lcarson.cols + j], ref.L, 1.0e-12);
      }
    }

    success *= close(model.y().getData()[Rcarson.offset + 1], model.y().getData()[Rcarson.offset + 2]);
    success *= close(model.y().getData()[Lcarson.offset + 1], model.y().getData()[Lcarson.offset + 2]);

    model.evaluateResidual();
    for (index_type i = 0; i < model.carson().size(); ++i)
    {
      success *= close(model.getResidual().getData()[eta.offset + i], 0.0, 1.0e-10);
    }

    success *= jacobianMatchesFiniteDifference(model, eta.offset, eta.offset, 1.0, 1.0e-7);
    success *= jacobianMatchesFiniteDifference(model, gammaR.offset, eta.offset, 1.0, 1.0e-7);
    success *= jacobianMatchesFiniteDifference(model, gammaR.offset, gammaR.offset, 1.0, 1.0e-7);
    success *= jacobianMatchesFiniteDifference(model, Rcarson.offset, gammaR.offset, 1.0, 1.0e-4);
    success *= jacobianMatchesFiniteDifference(model, Lcarson.offset, gammaI.offset, 1.0, 1.0e-4);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome series_impedance_model()
  {
    using GridKit::Testing::TestStatus;

    GridKit::EMT::Parameters::Overhead<scalar_type, index_type> model(makeData());
    model.updateTime(2.0 * pi * 50.0, 1.0);
    model.allocate();

    const auto rskin   = model.skinEffect().rskin();
    const auto lskin   = model.skinEffect().lskin();
    const auto Lgeo    = model.geometricInductance().Lgeo();
    const auto Rcarson = model.carson().Rcarson();
    const auto Lcarson = model.carson().Lcarson();
    const auto R       = model.seriesImpedance().R();
    const auto L       = model.seriesImpedance().L();

    TestStatus success  = true;
    success            *= (model.seriesImpedance().size() == 24);
    success            *= (R.rows == 2);
    success            *= (R.cols == 2);
    success            *= (L.rows == 2);
    success            *= (L.cols == 2);

    for (index_type i = 0; i < R.rows; ++i)
    {
      for (index_type j = 0; j < R.cols; ++j)
      {
        const index_type  k    = i * R.cols + j;
        const scalar_type Rref = ((i == j) ? model.y().getData()[rskin.offset + i] : 0.0)
                                 + model.y().getData()[Rcarson.offset + k];
        const scalar_type Lref = ((i == j) ? model.y().getData()[lskin.offset + i] : 0.0)
                                 + model.y().getData()[Lgeo.offset + k]
                                 + model.y().getData()[Lcarson.offset + k];
        success *= close(model.y().getData()[R.offset + k], Rref, 1.0e-12);
        success *= close(model.y().getData()[L.offset + k], Lref, 1.0e-12);
      }
    }

    model.evaluateResidual();
    const index_type series_base = R.offset - 4 * R.size();
    for (index_type i = 0; i < model.seriesImpedance().size(); ++i)
    {
      success *= close(model.getResidual().getData()[series_base + i], 0.0, 1.0e-10);
    }

    const index_type KK   = R.size();
    const index_type Rint = series_base;
    const index_type Lint = series_base + KK;
    const index_type Rext = series_base + 2 * KK;
    const index_type Lext = series_base + 3 * KK;

    success *= jacobianMatchesFiniteDifference(model, Rint, Rint, 1.0, 1.0e-7);
    success *= jacobianMatchesFiniteDifference(model, Rint, rskin.offset, 1.0, 1.0e-7);
    success *= jacobianMatchesFiniteDifference(model, Lint, lskin.offset, 1.0, 1.0e-7);
    success *= jacobianMatchesFiniteDifference(model, Rext, Rcarson.offset, 1.0, 1.0e-7);
    success *= jacobianMatchesFiniteDifference(model, Lext, Lgeo.offset, 1.0, 1.0e-7);
    success *= jacobianMatchesFiniteDifference(model, Lext, Lcarson.offset, 1.0, 1.0e-7);
    success *= jacobianMatchesFiniteDifference(model, R.offset, Rint, 1.0, 1.0e-7);
    success *= jacobianMatchesFiniteDifference(model, L.offset, Lext, 1.0, 1.0e-7);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome shunt_potential_scalar_model()
  {
    using GridKit::Testing::TestStatus;

    const scalar_type omega = 2.0 * pi * 50.0;
    auto              data  = makeScalarData();

    GridKit::EMT::Parameters::Conductor<scalar_type, index_type> conductor(data.conductor);
    conductor.initialize();
    GridKit::EMT::Parameters::Tower<scalar_type, index_type>          tower(data.tower, conductor);
    GridKit::EMT::Parameters::ShuntPotential<scalar_type, index_type> shunt(tower, conductor);

    tower.initialize();

    std::vector<scalar_type> y(static_cast<size_t>(shunt.size()), 0.0);
    std::vector<scalar_type> yp(static_cast<size_t>(shunt.size()), 0.0);
    std::vector<scalar_type> f(static_cast<size_t>(shunt.size()), 0.0);
    std::vector<index_type>  index(static_cast<size_t>(shunt.size()), 0);
    std::iota(index.begin(), index.end(), index_type{0});

    shunt.bind(y.data(), yp.data(), f.data(), index.data(), 0);
    shunt.setCoordinate(omega, 1.0);
    shunt.initialize();

    const auto        Gpot  = shunt.Gpot();
    const auto        Cpot  = shunt.Cpot();
    const scalar_type P_ref = std::log(2.0 * data.tower.height[0] / data.conductor.radius[0]) / (2.0 * pi * eps0);
    const scalar_type C_ref = 1.0 / P_ref;

    TestStatus success  = true;
    success            *= (shunt.size() == 3);
    success            *= (Gpot.offset == 1);
    success            *= (Cpot.offset == 2);
    success            *= close(y[0], P_ref);
    success            *= close(y[Gpot.offset], 0.0);
    success            *= close(y[Cpot.offset], C_ref);
    success            *= (y[Cpot.offset] > 0.0);

    shunt.evaluateResidual();
    for (const auto value : f)
    {
      success *= close(value, 0.0, 1.0e-9);
    }

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome shunt_potential_two_conductor_model()
  {
    using GridKit::Testing::TestStatus;

    GridKit::EMT::Parameters::Overhead<scalar_type, index_type> model(makeData());
    model.updateTime(2.0 * pi * 50.0, 1.0);
    model.allocate();

    const auto       Rcarson     = model.carson().Rcarson();
    const auto       Gpot        = model.shuntPotential().Gpot();
    const auto       Cpot        = model.shuntPotential().Cpot();
    const index_type carson_base = Rcarson.offset;

    TestStatus success  = true;
    success            *= (Gpot.rows == 2);
    success            *= (Gpot.cols == 2);
    success            *= (Cpot.rows == 2);
    success            *= (Cpot.cols == 2);
    for (index_type i = 0; i < Gpot.size(); ++i)
    {
      success *= close(model.y().getData()[Gpot.offset + i], 0.0);
    }
    success *= close(model.y().getData()[Gpot.offset + 1], model.y().getData()[Gpot.offset + 2]);
    success *= close(model.y().getData()[Cpot.offset + 1], model.y().getData()[Cpot.offset + 2]);

    model.evaluateResidual();
    for (index_type i = carson_base; i < model.size(); ++i)
    {
      success *= close(model.getResidual().getData()[i], 0.0, 1.0e-8);
    }

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome shunt_potential_three_conductor_passivity_model()
  {
    using GridKit::Testing::TestStatus;

    GridKit::EMT::Parameters::Overhead<scalar_type, index_type> model(makeThreeConductorData());
    model.updateTime(2.0 * pi * 1.0e8, 1.0);
    model.allocate();

    const auto Gpot = model.shuntPotential().Gpot();
    const auto Cpot = model.shuntPotential().Cpot();

    auto value = [&](const auto& signal, index_type i, index_type j)
    {
      return model.y().getData()[signal.offset + i * signal.cols + j];
    };

    std::array<scalar_type, 9> Cmat{};
    for (index_type i = 0; i < Cpot.rows; ++i)
    {
      for (index_type j = 0; j < Cpot.cols; ++j)
      {
        Cmat[i * Cpot.cols + j] = value(Cpot, i, j);
      }
    }

    const scalar_type det2 = Cmat[0] * Cmat[4] - Cmat[1] * Cmat[3];

    TestStatus success = true;
    for (index_type i = 0; i < Gpot.size(); ++i)
    {
      success *= close(model.y().getData()[Gpot.offset + i], 0.0);
    }
    for (index_type i = 0; i < Cpot.rows; ++i)
    {
      for (index_type j = 0; j < Cpot.cols; ++j)
      {
        success *= close(value(Cpot, i, j), value(Cpot, j, i), 1.0e-10);
      }
    }
    success *= (Cmat[0] > 0.0);
    success *= (det2 > 0.0);
    success *= (determinant3(Cmat) > 0.0);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome shunt_admittance_model()
  {
    using GridKit::Testing::TestStatus;

    GridKit::EMT::Parameters::Overhead<scalar_type, index_type> model(makeData());
    model.updateTime(2.0 * pi * 50.0, 1.0);
    model.allocate();

    const auto Gpot = model.shuntPotential().Gpot();
    const auto Cpot = model.shuntPotential().Cpot();
    const auto G    = model.shuntAdmittance().G();
    const auto C    = model.shuntAdmittance().C();

    TestStatus success  = true;
    success            *= (model.shuntAdmittance().size() == 16);
    success            *= (G.rows == 2);
    success            *= (G.cols == 2);
    success            *= (C.rows == 2);
    success            *= (C.cols == 2);

    for (index_type i = 0; i < G.size(); ++i)
    {
      success *= close(model.y().getData()[G.offset + i], model.y().getData()[Gpot.offset + i], 1.0e-12);
      success *= close(model.y().getData()[C.offset + i], model.y().getData()[Cpot.offset + i], 1.0e-12);
    }

    model.evaluateResidual();
    const index_type KK         = G.size();
    const index_type shunt_base = G.offset - 2 * KK;
    const index_type Gext       = shunt_base;
    const index_type Cext       = shunt_base + KK;
    for (index_type i = 0; i < KK; ++i)
    {
      success *= close(model.y().getData()[Gext + i], 0.0);
      success *= close(model.y().getData()[Cext + i], 0.0);
    }
    for (index_type i = 0; i < model.shuntAdmittance().size(); ++i)
    {
      success *= close(model.getResidual().getData()[shunt_base + i], 0.0, 1.0e-10);
    }

    success *= jacobianMatchesFiniteDifference(model, Gext, Gext, 1.0, 1.0e-7);
    success *= jacobianMatchesFiniteDifference(model, Cext, Cext, 1.0, 1.0e-7);
    success *= jacobianMatchesFiniteDifference(model, G.offset, Gpot.offset, 1.0, 1.0e-7);
    success *= jacobianMatchesFiniteDifference(model, G.offset, Gext, 1.0, 1.0e-7);
    success *= jacobianMatchesFiniteDifference(model, C.offset, Cpot.offset, 1.0, 1.0e-7);
    success *= jacobianMatchesFiniteDifference(model, C.offset, Cext, 1.0, 1.0e-7);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome yc_model()
  {
    using GridKit::Testing::TestStatus;

    const scalar_type                                           omega = 2.0 * pi * 50.0;
    auto                                                        data  = makeData();
    GridKit::EMT::Parameters::Overhead<scalar_type, index_type> model(data);
    model.updateTime(omega, 1.0);
    model.allocate();

    const auto Gc = model.yc().Gc();
    const auto Bc = model.yc().Bc();

    TestStatus success  = true;
    success            *= (model.yc().size() == 8);
    success            *= (Gc.rows == 2);
    success            *= (Gc.cols == 2);
    success            *= (Bc.rows == 2);
    success            *= (Bc.cols == 2);

    model.evaluateResidual();
    for (index_type i = 0; i < model.yc().size(); ++i)
    {
      success *= close(model.getResidual().getData()[Gc.offset + i], 0.0, 1.0e-8);
    }

    // Independent physical reference: the characteristic admittance
    // solves the sandwich equation Yc Z Yc = Y and is symmetric for the
    // symmetric per-unit-length inputs; the model's own residual cannot
    // certify either.
    {
      const auto R = model.seriesImpedance().R();
      const auto L = model.seriesImpedance().L();
      const auto G = model.shuntAdmittance().G();
      const auto C = model.shuntAdmittance().C();

      const auto* y = model.y().getData();

      Complex2x2 Z;
      Complex2x2 Y;
      Complex2x2 Yc;
      for (index_type e = 0; e < 4; ++e)
      {
        Z[e]  = {y[R.offset + e], omega * y[L.offset + e]};
        Y[e]  = {y[G.offset + e], omega * y[C.offset + e]};
        Yc[e] = {y[Gc.offset + e], y[Bc.offset + e]};
      }

      const auto sandwich = multiply2(multiply2(Yc, Z), Yc);
      const auto scale    = maxAbs2(Y);
      for (index_type e = 0; e < 4; ++e)
      {
        success *= (std::abs(sandwich[e] - Y[e]) <= 1.0e-10 * scale);
      }
      success *= (std::abs(Yc[1] - Yc[2]) <= 1.0e-10 * maxAbs2(Yc));
    }

    success *= jacobianMatchesFiniteDifference(model, Gc.offset + 1, Gc.offset + 0, 1.0, 1.0e-5);
    success *= jacobianMatchesFiniteDifference(model, Bc.offset + 1, Bc.offset + 3, 1.0, 1.0e-5);

    // Exact sparsity: the sandwich residual entry (0,0) reads only row 0
    // and column 0 of the admittance blocks, so a coupling to the (1,1)
    // entry must be structurally absent, not merely small.
    model.evaluateJacobian();
    success *= !csrHasEntry(model.getCsrJacobian(), Gc.offset + 0, Gc.offset + 3);
    success *= !csrHasEntry(model.getCsrJacobian(), Gc.offset + 0, Bc.offset + 3);

    GridKit::EMT::Parameters::Overhead<scalar_type, index_type> scalar_model(makeScalarData());
    scalar_model.updateTime(omega, 1.0);
    scalar_model.allocate();

    const auto R1  = scalar_model.seriesImpedance().R();
    const auto L1  = scalar_model.seriesImpedance().L();
    const auto G1  = scalar_model.shuntAdmittance().G();
    const auto C1  = scalar_model.shuntAdmittance().C();
    const auto Gc1 = scalar_model.yc().Gc();
    const auto Bc1 = scalar_model.yc().Bc();

    const std::complex<scalar_type> Z{scalar_model.y().getData()[R1.offset], omega * scalar_model.y().getData()[L1.offset]};
    const std::complex<scalar_type> Y{scalar_model.y().getData()[G1.offset], omega * scalar_model.y().getData()[C1.offset]};
    const auto                      Ycref = std::sqrt(Y / Z);

    success *= close(scalar_model.y().getData()[Gc1.offset], Ycref.real(), 1.0e-8);
    success *= close(scalar_model.y().getData()[Bc1.offset], Ycref.imag(), 1.0e-8);

    scalar_model.evaluateResidual();
    for (index_type i = 0; i < scalar_model.yc().size(); ++i)
    {
      success *= close(scalar_model.getResidual().getData()[Gc1.offset + i], 0.0, 1.0e-8);
    }

    success *= jacobianMatchesFiniteDifference(scalar_model, Gc1.offset, G1.offset, 1.0, 1.0e-5);
    success *= jacobianMatchesFiniteDifference(scalar_model, Bc1.offset, C1.offset, 1.0, 1.0e-5);
    success *= jacobianMatchesFiniteDifference(scalar_model, Gc1.offset, R1.offset, 1.0, 1.0e-5);
    success *= jacobianMatchesFiniteDifference(scalar_model, Gc1.offset, Gc1.offset, 1.0, 1.0e-5);
    success *= jacobianMatchesFiniteDifference(scalar_model, Bc1.offset, Bc1.offset, 1.0, 1.0e-5);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome zc_model()
  {
    using GridKit::Testing::TestStatus;

    const scalar_type                                           omega = 2.0 * pi * 50.0;
    auto                                                        data  = makeData();
    GridKit::EMT::Parameters::Overhead<scalar_type, index_type> model(data);
    model.updateTime(omega, 1.0);
    model.allocate();

    const auto Rc = model.zc().Rc();
    const auto Xc = model.zc().Xc();

    TestStatus success  = true;
    success            *= (model.zc().size() == 8);
    success            *= (Rc.rows == 2);
    success            *= (Rc.cols == 2);
    success            *= (Xc.rows == 2);
    success            *= (Xc.cols == 2);

    model.evaluateResidual();
    for (index_type i = 0; i < model.zc().size(); ++i)
    {
      success *= close(model.getResidual().getData()[Rc.offset + i], 0.0, 1.0e-8);
    }

    // Independent physical reference: the characteristic impedance
    // solves the sandwich equation Zc Y Zc = Z and is symmetric for the
    // symmetric per-unit-length inputs; the model's own residual cannot
    // certify either.
    {
      const auto R = model.seriesImpedance().R();
      const auto L = model.seriesImpedance().L();
      const auto G = model.shuntAdmittance().G();
      const auto C = model.shuntAdmittance().C();

      const auto* y = model.y().getData();

      Complex2x2 Z;
      Complex2x2 Y;
      Complex2x2 Zc;
      for (index_type e = 0; e < 4; ++e)
      {
        Z[e]  = {y[R.offset + e], omega * y[L.offset + e]};
        Y[e]  = {y[G.offset + e], omega * y[C.offset + e]};
        Zc[e] = {y[Rc.offset + e], y[Xc.offset + e]};
      }

      const auto sandwich = multiply2(multiply2(Zc, Y), Zc);
      const auto scale    = maxAbs2(Z);
      for (index_type e = 0; e < 4; ++e)
      {
        success *= (std::abs(sandwich[e] - Z[e]) <= 1.0e-10 * scale);
      }
      success *= (std::abs(Zc[1] - Zc[2]) <= 1.0e-10 * maxAbs2(Zc));
    }

    success *= jacobianMatchesFiniteDifference(model, Rc.offset + 1, Rc.offset + 0, 1.0, 1.0e-5);
    success *= jacobianMatchesFiniteDifference(model, Xc.offset + 1, Xc.offset + 3, 1.0, 1.0e-5);

    // Exact sparsity: the sandwich residual entry (0,0) reads only row 0
    // and column 0 of the impedance blocks, so a coupling to the (1,1)
    // entry must be structurally absent, not merely small.
    model.evaluateJacobian();
    success *= !csrHasEntry(model.getCsrJacobian(), Rc.offset + 0, Rc.offset + 3);
    success *= !csrHasEntry(model.getCsrJacobian(), Rc.offset + 0, Xc.offset + 3);

    GridKit::EMT::Parameters::Overhead<scalar_type, index_type> scalar_model(makeScalarData());
    scalar_model.updateTime(omega, 1.0);
    scalar_model.allocate();

    const auto R1  = scalar_model.seriesImpedance().R();
    const auto L1  = scalar_model.seriesImpedance().L();
    const auto G1  = scalar_model.shuntAdmittance().G();
    const auto C1  = scalar_model.shuntAdmittance().C();
    const auto Rc1 = scalar_model.zc().Rc();
    const auto Xc1 = scalar_model.zc().Xc();

    const std::complex<scalar_type> Z{scalar_model.y().getData()[R1.offset], omega * scalar_model.y().getData()[L1.offset]};
    const std::complex<scalar_type> Y{scalar_model.y().getData()[G1.offset], omega * scalar_model.y().getData()[C1.offset]};
    const auto                      Zcref = std::sqrt(Z / Y);

    success *= close(scalar_model.y().getData()[Rc1.offset], Zcref.real(), 1.0e-8);
    success *= close(scalar_model.y().getData()[Xc1.offset], Zcref.imag(), 1.0e-8);

    scalar_model.evaluateResidual();
    for (index_type i = 0; i < scalar_model.zc().size(); ++i)
    {
      success *= close(scalar_model.getResidual().getData()[Rc1.offset + i], 0.0, 1.0e-8);
    }

    success *= jacobianMatchesFiniteDifference(scalar_model, Rc1.offset, R1.offset, 1.0, 1.0e-5);
    success *= jacobianMatchesFiniteDifference(scalar_model, Xc1.offset, L1.offset, 1.0, 1.0e-5);
    success *= jacobianMatchesFiniteDifference(scalar_model, Rc1.offset, G1.offset, 1.0, 1.0e-5);
    success *= jacobianMatchesFiniteDifference(scalar_model, Rc1.offset, Rc1.offset, 1.0, 1.0e-5);
    success *= jacobianMatchesFiniteDifference(scalar_model, Xc1.offset, Xc1.offset, 1.0, 1.0e-5);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome gamma_model()
  {
    using GridKit::Testing::TestStatus;

    const scalar_type                                           omega = 2.0 * pi * 50.0;
    GridKit::EMT::Parameters::Overhead<scalar_type, index_type> model(makeData());
    model.updateTime(omega, 1.0);
    model.allocate();

    const auto alpha = model.gamma().alpha();
    const auto beta  = model.gamma().beta();

    TestStatus success  = true;
    success            *= (model.gamma().size() == 8);
    success            *= (alpha.rows == 2 && alpha.cols == 2);
    success            *= (beta.rows == 2 && beta.cols == 2);

    model.evaluateResidual();
    for (index_type i = 0; i < model.gamma().size(); ++i)
    {
      success *= close(model.getResidual().getData()[alpha.offset + i], 0.0, 1.0e-8);
    }

    GridKit::EMT::Parameters::Overhead<scalar_type, index_type> scalar_model(makeScalarData());
    scalar_model.updateTime(omega, 1.0);
    scalar_model.allocate();

    const auto R1 = scalar_model.seriesImpedance().R();
    const auto L1 = scalar_model.seriesImpedance().L();
    const auto G1 = scalar_model.shuntAdmittance().G();
    const auto C1 = scalar_model.shuntAdmittance().C();
    const auto ar = scalar_model.gamma().alpha();
    const auto bi = scalar_model.gamma().beta();

    const std::complex<scalar_type> Z{scalar_model.y().getData()[R1.offset], omega * scalar_model.y().getData()[L1.offset]};
    const std::complex<scalar_type> Y{scalar_model.y().getData()[G1.offset], omega * scalar_model.y().getData()[C1.offset]};
    const auto                      gamma_ref = std::sqrt(Z * Y);

    success *= (scalar_model.gamma().size() == 2);
    success *= close(scalar_model.y().getData()[ar.offset], gamma_ref.real(), 1.0e-8);
    success *= close(scalar_model.y().getData()[bi.offset], gamma_ref.imag(), 1.0e-8);

    scalar_model.evaluateResidual();
    for (index_type i = 0; i < scalar_model.gamma().size(); ++i)
    {
      success *= close(scalar_model.getResidual().getData()[ar.offset + i], 0.0, 1.0e-8);
    }

    success *= jacobianMatchesFiniteDifference(scalar_model, ar.offset, ar.offset, 1.0, 1.0e-5);
    success *= jacobianMatchesFiniteDifference(scalar_model, bi.offset, bi.offset, 1.0, 1.0e-5);
    success *= jacobianMatchesFiniteDifference(scalar_model, ar.offset, R1.offset, 1.0, 1.0e-5);
    success *= jacobianMatchesFiniteDifference(scalar_model, bi.offset, C1.offset, 1.0, 1.0e-5);

    return success.report(__func__);
  }

  /// Modal quantities are observations refreshed at emission, not DAE
  /// states: after printMonitoredVariables() the decomposition must be
  /// biorthogonal, satisfy the eigen-relation against the integrated
  /// gamma states, and carry the closed-form modal responses.
  GridKit::Testing::TestOutcome overhead_modal_observation()
  {
    using GridKit::Testing::TestStatus;
    using MonitorVariable =
        GridKit::EMT::Parameters::OverheadData<scalar_type, index_type>::MonitorableVariables;
    using ComplexT = std::complex<scalar_type>;

    const scalar_type omega  = 2.0 * pi * 50.0;
    auto              data   = makeData();
    data.monitored_variables = {MonitorVariable::Tv,
                                MonitorVariable::Ti,
                                MonitorVariable::Alpha,
                                MonitorVariable::Beta,
                                MonitorVariable::Tau,
                                MonitorVariable::H};

    GridKit::EMT::Parameters::Overhead<scalar_type, index_type> model(data);
    model.updateTime(omega, 1.0);
    model.allocate();

    TestStatus success  = true;
    success            *= (model.modes() != nullptr);
    if (model.modes() == nullptr)
    {
      return success.report(__func__);
    }

    // Two emissions: the first fixes the canonical frame, the second
    // exercises the identity-tracked path.
    model.printMonitoredVariables();
    model.printMonitoredVariables();
    model.stopMonitor();

    const auto&      modes   = *model.modes();
    const auto*      values  = modes.data();
    const auto       ell     = model.path().length();
    const index_type k       = modes.modeCount();
    success                 *= (k == 2);

    // Biorthogonality is exact by construction.
    for (index_type m = 0; m < k; ++m)
    {
      for (index_type n = 0; n < k; ++n)
      {
        ComplexT pairing{0.0, 0.0};
        for (index_type i = 0; i < k; ++i)
        {
          const auto     entry = i * k;
          const ComplexT tv{values[modes.tvReal().offset + entry + n],
                            values[modes.tvImag().offset + entry + n]};
          const ComplexT ti{values[modes.tiReal().offset + entry + m],
                            values[modes.tiImag().offset + entry + m]};
          pairing += std::conj(ti) * tv;
        }
        success *= close(pairing.real(), (m == n) ? 1.0 : 0.0, 1.0e-12);
        success *= close(pairing.imag(), 0.0, 1.0e-12);
      }
    }

    // The eigen-relation certifies the observation against the
    // integrated gamma states.
    const auto alpha = model.gamma().alpha();
    const auto beta  = model.gamma().beta();
    for (index_type m = 0; m < k; ++m)
    {
      const ComplexT lambda{values[modes.alphaM().offset + m],
                            values[modes.betaM().offset + m]};
      for (index_type i = 0; i < k; ++i)
      {
        ComplexT product{0.0, 0.0};
        for (index_type j = 0; j < k; ++j)
        {
          const ComplexT gamma_entry{
              model.y().getData()[alpha.offset + i * k + j],
              model.y().getData()[beta.offset + i * k + j]};
          const ComplexT tv{values[modes.tvReal().offset + j * k + m],
                            values[modes.tvImag().offset + j * k + m]};
          product += gamma_entry * tv;
        }
        const ComplexT scaled =
            lambda
            * ComplexT{values[modes.tvReal().offset + i * k + m],
                       values[modes.tvImag().offset + i * k + m]};
        success *= close(product.real(), scaled.real(), 1.0e-8);
        success *= close(product.imag(), scaled.imag(), 1.0e-8);
      }

      // Closed-form modal responses.
      const ComplexT factor  = std::exp(-ell * lambda);
      success               *= close(values[modes.hReal().offset + m], factor.real(), 1.0e-10);
      success               *= close(values[modes.hImag().offset + m], factor.imag(), 1.0e-10);
      success               *= close(values[modes.tau().offset + m],
                       ell * lambda.imag() / omega,
                       1.0e-10);
    }

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome overhead_gis_path_length()
  {
    using GridKit::Testing::TestStatus;

    const scalar_type omega = 2.0 * pi * 50.0;
    auto              data  = makeData();
    data.path.length.reset();
    data.path.path = {{0.0, 0.0}, {0.0, 1.0}};

    GridKit::EMT::Parameters::Overhead<scalar_type, index_type> model(data);
    model.updateTime(omega, 1.0);
    model.allocate();

    const scalar_type ell_ref =
        GridKit::EMT::Parameters::Path<scalar_type, index_type>::earthRadius() * pi / 180.0
        * model.tower().saggedToSpanRatio();
    const scalar_type ell = model.path().length();

    TestStatus success  = true;
    success            *= close(ell, ell_ref, 1.0e-12);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome overhead_layout()
  {
    using GridKit::Testing::TestStatus;

    const scalar_type omega = 2.0 * pi * 50.0;

    GridKit::EMT::Parameters::Overhead<scalar_type, index_type> model(makeData());
    model.updateTime(omega, 1.0);
    model.allocate();

    const auto eta     = model.carson().eta();
    const auto gammaR  = model.carson().gammaR();
    const auto gammaI  = model.carson().gammaI();
    const auto Lgeo    = model.geometricInductance().Lgeo();
    const auto rskin   = model.skinEffect().rskin();
    const auto lskin   = model.skinEffect().lskin();
    const auto Rcarson = model.carson().Rcarson();
    const auto Lcarson = model.carson().Lcarson();
    const auto R       = model.seriesImpedance().R();
    const auto L       = model.seriesImpedance().L();
    const auto Gpot    = model.shuntPotential().Gpot();
    const auto Cpot    = model.shuntPotential().Cpot();
    const auto G       = model.shuntAdmittance().G();
    const auto C       = model.shuntAdmittance().C();
    const auto alpha   = model.gamma().alpha();
    const auto beta    = model.gamma().beta();
    const auto Gc      = model.yc().Gc();
    const auto Bc      = model.yc().Bc();
    const auto Rc      = model.zc().Rc();
    const auto Xc      = model.zc().Xc();

    TestStatus       success = true;
    const index_type skin_size =
        2 * GridKit::EMT::Parameters::SkinEffect<scalar_type, index_type>::rootCount()
        + 5 * model.conductor().conductorCount();
    success *= (model.tower().size() == 0);
    success *= (model.conductor().size() == 0);
    success *= (model.path().size() == 0);
    success *= close(model.path().length(), 100000.0);
    success *= close(model.tower().span(), 300.0);
    success *= (model.conductor().conductorCount() == 2);
    success *= close(model.conductor().weight(0), 11.0);
    success *= (model.conductor().phase(0) == "a");
    success *= (model.geometricInductance().size() == 4);
    success *= (model.skinEffect().size() == skin_size);
    success *= (model.carson().size() == 11);
    success *= (model.seriesImpedance().size() == 24);
    success *= (model.shuntPotential().size() == 12);
    success *= (model.shuntAdmittance().size() == 16);
    success *= (model.gamma().size() == 8);
    success *= (model.yc().size() == 8);
    success *= (model.zc().size() == 8);
    success *= (Lgeo.offset == 0);
    success *= (rskin.offset == Lgeo.offset + Lgeo.size() + model.skinEffect().size() - rskin.size() - lskin.size());
    success *= (lskin.offset == rskin.offset + rskin.size());
    success *= (eta.offset == lskin.offset + lskin.size());
    success *= (gammaR.offset == eta.offset + eta.size());
    success *= (gammaI.offset == gammaR.offset + gammaR.size());
    success *= (Rcarson.offset == gammaI.offset + gammaI.size());
    success *= (Lcarson.offset == Rcarson.offset + Rcarson.size());
    success *= (R.offset == Lcarson.offset + Lcarson.size() + 4 * R.size());
    success *= (L.offset == R.offset + R.size());
    success *= (Gpot.offset == L.offset + L.size() + Gpot.size());
    success *= (Cpot.offset == Gpot.offset + Gpot.size());
    success *= (G.offset == Cpot.offset + Cpot.size() + 2 * G.size());
    success *= (C.offset == G.offset + G.size());
    success *= (alpha.offset == C.offset + C.size());
    success *= (beta.offset == alpha.offset + alpha.size());
    success *= (Gc.offset == beta.offset + beta.size());
    success *= (Bc.offset == Gc.offset + Gc.size());
    success *= (Rc.offset == Bc.offset + Bc.size());
    success *= (Xc.offset == Rc.offset + Rc.size());
    success *= (model.size() == Xc.offset + Xc.size());

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome overhead_monitor_outputs()
  {
    using GridKit::Testing::TestStatus;
    using MonitorVariable =
        GridKit::EMT::Parameters::OverheadData<scalar_type, index_type>::MonitorableVariables;

    const std::string path =
        (std::filesystem::temp_directory_path() / "gridkit_overhead_monitor_test.csv").string();
    std::remove(path.c_str());

    auto data                = makeData();
    data.monitored_variables = {MonitorVariable::R,
                                MonitorVariable::L,
                                MonitorVariable::G,
                                MonitorVariable::C,
                                MonitorVariable::Tv,
                                MonitorVariable::Ti,
                                MonitorVariable::Alpha,
                                MonitorVariable::Beta,
                                MonitorVariable::Tau,
                                MonitorVariable::H,
                                MonitorVariable::Yc,
                                MonitorVariable::Zc};
    data.monitor_sink        = {{GridKit::Model::VariableMonitorFormat::CSV, path}};

    GridKit::EMT::Parameters::Overhead<scalar_type, index_type> model(data);
    model.updateTime(2.0 * pi * 50.0, 1.0);
    model.allocate();
    model.printMonitoredVariables();
    model.stopMonitor();

    std::ifstream stream(path);
    std::string   header;
    std::getline(stream, header);
    std::remove(path.c_str());

    TestStatus success = true;
    for (const auto* name : {"Overhead_R_0_0",
                             "Overhead_L_0_0",
                             "Overhead_G_0_0",
                             "Overhead_C_0_0",
                             "Overhead_Tv_real_0_0",
                             "Overhead_Tv_imag_0_0",
                             "Overhead_Ti_real_0_0",
                             "Overhead_Ti_imag_0_0",
                             "Overhead_Alpha_0",
                             "Overhead_Beta_0",
                             "Overhead_Tau_0",
                             "Overhead_H_real_0",
                             "Overhead_H_imag_0",
                             "Overhead_Yc_real_0_0",
                             "Overhead_Yc_imag_0_0",
                             "Overhead_Zc_real_0_0",
                             "Overhead_Zc_imag_0_0"})
    {
      success *= (header.find(name) != std::string::npos);
    }

    success *= (header.find("Rprime") == std::string::npos);
    success *= (header.find("Lprime") == std::string::npos);
    success *= (header.find("Overhead_X_0_0") == std::string::npos);
    success *= (header.find("Overhead_B_0_0") == std::string::npos);
    success *= (header.find("Bprime") == std::string::npos);
    success *= (header.find("ARight") == std::string::npos);
    success *= (header.find("BRight") == std::string::npos);
    success *= (header.find("ALeft") == std::string::npos);
    success *= (header.find("BLeft") == std::string::npos);
    success *= (header.find("Overhead_Alpha_0_0") == std::string::npos);
    success *= (header.find("Overhead_Beta_0_0") == std::string::npos);
    success *= (header.find("Overhead_H_real_0_0") == std::string::npos);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome reduction_map_model()
  {
    using GridKit::Testing::TestStatus;
    using MapT = GridKit::EMT::Parameters::ReductionMap<scalar_type, index_type>;

    TestStatus success = true;

    // Distinct phases and no grounding: the identity, so the aggregate
    // skips the reduction stage entirely.
    const auto plain  = MapT::fromConductors(makeData().conductor);
    success          *= plain.identity();
    success          *= (plain.terminals == 2);
    success          *= (plain.terminal[0] == 0 && plain.terminal[1] == 1);

    const auto bundled  = MapT::fromConductors(makeBundledData().conductor);
    success            *= !bundled.identity();
    success            *= (bundled.conductors == 5);
    success            *= (bundled.terminals == 2);
    success            *= (bundled.terminal[0] == 0 && bundled.terminal[1] == 0);
    success            *= (bundled.terminal[2] == 1 && bundled.terminal[3] == 1);
    success            *= (bundled.terminal[4] == MapT::kGrounded);
    success            *= (bundled.members[0] == std::vector<index_type>{0, 1});
    success            *= (bundled.members[1] == std::vector<index_type>{2, 3});

    // The same phase on different circuits is a different terminal.
    auto circuits     = makeData().conductor;
    circuits.phase    = {"a", "a"};
    circuits.circuit  = {1, 2};
    const auto split  = MapT::fromConductors(circuits);
    success          *= (split.terminals == 2);
    success          *= split.identity();

    // A line with every conductor grounded has no terminals.
    auto grounded  = makeData().conductor;
    grounded.phase = {"g", "g"};
    bool thrown    = false;
    try
    {
      MapT::fromConductors(grounded);
    }
    catch (const std::exception&)
    {
      thrown = true;
    }
    success *= thrown;

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome series_reduction_model()
  {
    using GridKit::Testing::TestStatus;
    using ComplexT = std::complex<scalar_type>;
    using MatrixXc = Eigen::MatrixXcd;

    const scalar_type                                           omega = 2.0 * pi * 50.0;
    auto                                                        data  = makeBundledData();
    GridKit::EMT::Parameters::Overhead<scalar_type, index_type> model(data);
    model.updateTime(omega, 1.0);
    model.allocate();

    TestStatus success  = true;
    success            *= (model.seriesReduction() != nullptr);
    if (model.seriesReduction() == nullptr)
    {
      return success.report(__func__);
    }

    const auto& reduction  = *model.seriesReduction();
    const auto  Drs        = reduction.Dr();
    const auto  Rr         = reduction.R();
    const auto  Lr         = reduction.L();
    success               *= (Rr.rows == 2 && Rr.cols == 2);

    // Exact initialization: every reduction residual row vanishes.
    const auto* f = model.getResidual().getData();
    for (index_type i = 0; i < 2 * 10 + 2 * 4; ++i)
    {
      success *= close(f[Drs.offset + i], 0.0, 1.0e-10);
    }

    // Independent dense reference through the explicit inverses the
    // element never forms: Kron-eliminate the grounded conductor, then
    // Zred = (E^T Zk^-1 E)^-1.
    const auto  R = model.seriesImpedance().R();
    const auto  L = model.seriesImpedance().L();
    const auto* y = model.y().getData();

    MatrixXc Z(5, 5);
    for (index_type i = 0; i < 5; ++i)
    {
      for (index_type j = 0; j < 5; ++j)
      {
        Z(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) =
            ComplexT{y[R.offset + i * 5 + j], omega * y[L.offset + i * 5 + j]};
      }
    }
    const MatrixXc Zk = Z.topLeftCorner(4, 4)
                        - Z.topRightCorner(4, 1)
                              * Z.bottomRightCorner(1, 1).inverse()
                              * Z.bottomLeftCorner(1, 4);
    MatrixXc incidence = MatrixXc::Zero(4, 2);
    incidence(0, 0) = incidence(1, 0) = ComplexT{1.0, 0.0};
    incidence(2, 1) = incidence(3, 1) = ComplexT{1.0, 0.0};
    const MatrixXc reference =
        (incidence.transpose() * Zk.inverse() * incidence).inverse();

    const scalar_type scale = reference.cwiseAbs().maxCoeff();
    for (index_type t = 0; t < 2; ++t)
    {
      for (index_type q = 0; q < 2; ++q)
      {
        const ComplexT value{y[Rr.offset + t * 2 + q],
                             omega * y[Lr.offset + t * 2 + q]};
        success *= (std::abs(value
                             - reference(static_cast<Eigen::Index>(t),
                                         static_cast<Eigen::Index>(q)))
                    <= 1.0e-10 * scale);
      }
    }

    // The current-sum rows have exact unit couplings to D, and the
    // distribution rows couple to Zred with exact -1.
    success *= jacobianMatchesFiniteDifference(model, Rr.offset, Drs.offset, 1.0, 1.0e-7);
    success *= jacobianMatchesFiniteDifference(model, Drs.offset, Rr.offset, 1.0, 1.0e-7);

    // Conductor 0 belongs to terminal 0, so its distribution rows must
    // be structurally blind to the terminal-1 impedance entries.
    model.evaluateJacobian();
    success *= !csrHasEntry(model.getCsrJacobian(), Drs.offset, Rr.offset + 3);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome shunt_reduction_model()
  {
    using GridKit::Testing::TestStatus;

    const scalar_type                                           omega = 2.0 * pi * 50.0;
    auto                                                        data  = makeBundledData();
    GridKit::EMT::Parameters::Overhead<scalar_type, index_type> model(data);
    model.updateTime(omega, 1.0);
    model.allocate();

    TestStatus success  = true;
    success            *= (model.shuntReduction() != nullptr);
    if (model.shuntReduction() == nullptr)
    {
      return success.report(__func__);
    }

    const auto& reduction  = *model.shuntReduction();
    const auto  Gr         = reduction.G();
    const auto  Cr         = reduction.C();
    success               *= (Gr.rows == 2 && Gr.cols == 2);

    const auto* f = model.getResidual().getData();
    for (index_type i = 0; i < 2 * 4; ++i)
    {
      success *= close(f[Gr.offset + i], 0.0, 1.0e-10);
    }

    // Independent congruence reference: Yred = E^T Y E is a plain sum
    // over bundle members.
    const auto                      G          = model.shuntAdmittance().G();
    const auto                      C          = model.shuntAdmittance().C();
    const auto*                     y          = model.y().getData();
    const std::array<index_type, 2> members[2] = {{0, 1}, {2, 3}};

    for (index_type t = 0; t < 2; ++t)
    {
      for (index_type q = 0; q < 2; ++q)
      {
        scalar_type sum_g = 0.0;
        scalar_type sum_c = 0.0;
        for (const index_type a : members[t])
        {
          for (const index_type b : members[q])
          {
            sum_g += y[G.offset + a * 5 + b];
            sum_c += y[C.offset + a * 5 + b];
          }
        }
        success *= close(y[Gr.offset + t * 2 + q], sum_g, 1.0e-12);
        success *= close(y[Cr.offset + t * 2 + q], sum_c, 1.0e-12);
      }
    }

    success *= jacobianMatchesFiniteDifference(model, Gr.offset, Gr.offset, 1.0, 1.0e-7);
    success *= jacobianMatchesFiniteDifference(model, Gr.offset, G.offset, 1.0, 1.0e-7);

    // The conductance rows never touch the capacitance states.
    model.evaluateJacobian();
    success *= !csrHasEntry(model.getCsrJacobian(), Gr.offset, Cr.offset);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome overhead_reduced_chain()
  {
    using GridKit::Testing::TestStatus;

    const scalar_type                                           omega = 2.0 * pi * 50.0;
    auto                                                        data  = makeBundledData();
    GridKit::EMT::Parameters::Overhead<scalar_type, index_type> model(data);
    model.updateTime(omega, 1.0);
    model.allocate();

    // The whole modal stage runs at the terminal count.
    const auto Gc = model.yc().Gc();
    const auto Bc = model.yc().Bc();

    TestStatus success  = true;
    success            *= (Gc.rows == 2 && Gc.cols == 2);
    success            *= (model.gamma().alpha().rows == 2);

    if (model.seriesReduction() == nullptr || model.shuntReduction() == nullptr)
    {
      success *= false;
      return success.report(__func__);
    }

    // The reduced characteristic admittance satisfies the sandwich
    // equation against the reduced per-unit-length pair.
    const auto  Rr = model.seriesReduction()->R();
    const auto  Lr = model.seriesReduction()->L();
    const auto  Gr = model.shuntReduction()->G();
    const auto  Cr = model.shuntReduction()->C();
    const auto* y  = model.y().getData();

    Complex2x2 Z;
    Complex2x2 Y;
    Complex2x2 Yc;
    for (index_type e = 0; e < 4; ++e)
    {
      Z[e]  = {y[Rr.offset + e], omega * y[Lr.offset + e]};
      Y[e]  = {y[Gr.offset + e], omega * y[Cr.offset + e]};
      Yc[e] = {y[Gc.offset + e], y[Bc.offset + e]};
    }

    const auto sandwich = multiply2(multiply2(Yc, Z), Yc);
    const auto scale    = maxAbs2(Y);
    for (index_type e = 0; e < 4; ++e)
    {
      success *= (std::abs(sandwich[e] - Y[e]) <= 1.0e-10 * scale);
    }
    success *= (std::abs(Yc[1] - Yc[2]) <= 1.0e-10 * maxAbs2(Yc));

    return success.report(__func__);
  }
} // namespace

int main()
{
  GridKit::Testing::TestingResults result;
  result += tower_model();
  result += geometric_inductance_model();
  result += carson_model();
  result += series_impedance_model();
  result += shunt_potential_scalar_model();
  result += shunt_potential_two_conductor_model();
  result += shunt_potential_three_conductor_passivity_model();
  result += shunt_admittance_model();
  result += yc_model();
  result += zc_model();
  result += gamma_model();
  result += overhead_modal_observation();
  result += overhead_gis_path_length();
  result += overhead_layout();
  result += overhead_monitor_outputs();
  result += reduction_map_model();
  result += series_reduction_model();
  result += shunt_reduction_model();
  result += overhead_reduced_chain();
  return result.summary();
}
