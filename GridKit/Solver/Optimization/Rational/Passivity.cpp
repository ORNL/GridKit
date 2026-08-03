/**
 * @file Passivity.cpp
 *
 * @brief Passivity assessment of fitted rational models.
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#include "Passivity.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace
{
  constexpr size_t SCREEN_POINTS    = 512;
  constexpr size_t BISECTION_STEPS  = 60;
  constexpr double SCREEN_RANGE_PAD = 10.0;
} // namespace

namespace GridKit
{
  namespace Optimization
  {
    namespace Detail
    {
      /**
       * @brief Smallest eigenvalue of the Hermitian part of Y(j omega).
       */
      template <typename scalar_type, typename index_type>
      int conductanceFloor(
          const RationalModel<scalar_type, index_type>&       model,
          typename GridKit::ScalarTraits<scalar_type>::RealT  omega,
          typename GridKit::ScalarTraits<scalar_type>::RealT& floor)
      {
        using RealT    = typename GridKit::ScalarTraits<scalar_type>::RealT;
        using ComplexT = std::complex<RealT>;
        using MatrixXc =
            Eigen::Matrix<ComplexT, Eigen::Dynamic, Eigen::Dynamic>;

        const auto ports = static_cast<Eigen::Index>(model.rows);

        MatrixXc admittance(ports, ports);
        for (Eigen::Index row = 0; row < ports; ++row)
        {
          for (Eigen::Index col = 0; col < ports; ++col)
          {
            admittance(row, col) =
                model.evaluate(omega,
                               static_cast<index_type>(row),
                               static_cast<index_type>(col));
          }
        }

        const MatrixXc hermitian =
            RealT{0.5} * (admittance + admittance.adjoint());
        const Eigen::SelfAdjointEigenSolver<MatrixXc> eigensolver(
            hermitian);
        if (eigensolver.info() != Eigen::Success)
        {
          return -2;
        }
        floor = eigensolver.eigenvalues()(0);
        return 0;
      }
    } // namespace Detail

    /**
     * @brief Assess passivity of a fitted admittance model.
     *
     * An unstable model is never passive, so stability gates the screen.
     * The conductance floor -- the smallest eigenvalue of the Hermitian
     * part of Y(j omega) -- is then sampled at DC, where a fitted
     * admittance classically violates first, on a logarithmic grid
     * spanning the union of the fitted band and the pole magnitudes, and
     * at the asymptotic limit through the constant term. Sign changes
     * inside the grid are refined by bisection; a violation still open at
     * the grid end extends to infinity when the constant term is itself
     * indefinite. The exact Hamiltonian certification arrives with the
     * state-space realization; between grid points this screen remains
     * resolution limited by design.
     */
    template <typename scalar_type, typename index_type>
    int assessPassivity(
        const RationalModel<scalar_type, index_type>&      model,
        PassivityReport<scalar_type, index_type>&          report,
        typename GridKit::ScalarTraits<scalar_type>::RealT omega_low,
        typename GridKit::ScalarTraits<scalar_type>::RealT omega_high)
    {
      using RealT = typename GridKit::ScalarTraits<scalar_type>::RealT;

      report.passive = false;
      report.stable  = false;
      report.violations.clear();

      if (model.rows != model.cols || model.rows < 1)
      {
        return -1;
      }
      if (!(omega_low > RealT{0}) || !(omega_high > omega_low))
      {
        return -1;
      }

      report.stable = model.isStable(RealT{0});
      if (!report.stable)
      {
        return 0;
      }

      RealT low  = omega_low;
      RealT high = omega_high;
      for (const auto& pole : model.poles)
      {
        const RealT scale = std::abs(pole);
        if (scale > RealT{0})
        {
          low  = std::min(low, scale);
          high = std::max(high, scale);
        }
      }
      low  = low / RealT{SCREEN_RANGE_PAD};
      high = high * RealT{SCREEN_RANGE_PAD};

      RealT dc_floor = RealT{0};
      if (Detail::conductanceFloor(model, RealT{0}, dc_floor) != 0)
      {
        return -2;
      }

      std::vector<RealT> grid(SCREEN_POINTS, RealT{0});
      std::vector<RealT> floors(SCREEN_POINTS, RealT{0});
      for (size_t k = 0; k < SCREEN_POINTS; ++k)
      {
        const RealT fraction = static_cast<RealT>(k) / static_cast<RealT>(SCREEN_POINTS - 1);
        grid[k]              = std::exp(std::log(low) + fraction * (std::log(high) - std::log(low)));
        if (Detail::conductanceFloor(model, grid[k], floors[k]) != 0)
        {
          return -2;
        }
      }

      // Refine a sign change of the conductance floor by bisection; the
      // arithmetic midpoint supports a zero bracket edge, so the DC
      // sample joins the state machine like any grid point. An
      // eigensolver failure propagates instead of truncating the search.
      const auto crossing = [&model](RealT inside, RealT outside, RealT& edge) -> int
      {
        for (size_t step = 0; step < BISECTION_STEPS; ++step)
        {
          const RealT middle = (inside + outside) / RealT{2};
          RealT       floor  = RealT{0};
          if (Detail::conductanceFloor(model, middle, floor) != 0)
          {
            return -2;
          }
          if (floor < RealT{0})
          {
            inside = middle;
          }
          else
          {
            outside = middle;
          }
        }
        edge = (inside + outside) / RealT{2};
        return 0;
      };

      using ReportT   = PassivityReport<scalar_type, index_type>;
      using Violation = typename ReportT::Violation;

      bool      violating = dc_floor < RealT{0};
      Violation band{RealT{0}, RealT{0}};

      for (size_t k = 0; k < SCREEN_POINTS; ++k)
      {
        const RealT previous = k == 0 ? RealT{0} : grid[k - 1];
        const bool  now      = floors[k] < RealT{0};
        if (now && !violating)
        {
          if (crossing(grid[k], previous, band.omega_start) != 0)
          {
            return -2;
          }
          violating = true;
        }
        else if (!now && violating)
        {
          if (crossing(previous, grid[k], band.omega_end) != 0)
          {
            return -2;
          }
          report.violations.push_back(band);
          violating = false;
        }
      }

      // Asymptotic limit: Y(j omega) tends to the constant term, whose
      // symmetric-part floor decides passivity beyond any finite grid.
      const RealT infinity         = std::numeric_limits<RealT>::infinity();
      RealT       asymptotic_floor = RealT{0};
      if (!model.d.empty())
      {
        using MatrixXr = Eigen::Matrix<RealT, Eigen::Dynamic, Eigen::Dynamic>;

        const auto ports = static_cast<Eigen::Index>(model.rows);
        MatrixXr   constant(ports, ports);
        for (Eigen::Index row = 0; row < ports; ++row)
        {
          for (Eigen::Index col = 0; col < ports; ++col)
          {
            constant(row, col) = model.d[static_cast<size_t>(row * ports + col)];
          }
        }
        const MatrixXr                                symmetric = RealT{0.5} * (constant + constant.transpose());
        const Eigen::SelfAdjointEigenSolver<MatrixXr> eigensolver(symmetric);
        if (eigensolver.info() != Eigen::Success)
        {
          return -2;
        }
        asymptotic_floor = eigensolver.eigenvalues()(0);
      }

      if (violating)
      {
        band.omega_end =
            asymptotic_floor < RealT{0} ? infinity : grid[SCREEN_POINTS - 1];
        report.violations.push_back(band);
      }
      else if (asymptotic_floor < RealT{0})
      {
        report.violations.push_back(
            Violation{grid[SCREEN_POINTS - 1], infinity});
      }

      report.passive = report.violations.empty();
      return 0;
    }

    template int
    assessPassivity<double, long int>(const RationalModel<double, long int>&,
                                      PassivityReport<double, long int>&,
                                      double,
                                      double);
    template int
    assessPassivity<double, size_t>(const RationalModel<double, size_t>&,
                                    PassivityReport<double, size_t>&,
                                    double,
                                    double);
  } // namespace Optimization
} // namespace GridKit
