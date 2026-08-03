```cpp
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

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace
{
  constexpr size_t SCREEN_POINTS     = 512;
  constexpr size_t BISECTION_STEPS   = 60;
  constexpr double SCREEN_RANGE_PAD  = 10.0;
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
          const RationalModel<scalar_type, index_type>& model,
          typename GridKit::ScalarTraits<scalar_type>::RealT omega,
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
     * Screens the smallest eigenvalue of the Hermitian part of Y(j omega)
     * on a logarithmic grid spanning the pole frequencies, then refines
     * every sign change by bisection to locate violation band edges. The
     * exact Hamiltonian certification arrives with the state-space
     * realization; this screen is grid-resolution limited by design.
     */
    template <typename scalar_type, typename index_type>
    int assessPassivity(const RationalModel<scalar_type, index_type>& model,
                        PassivityReport<scalar_type, index_type>& report)
    {
      using RealT = typename GridKit::ScalarTraits<scalar_type>::RealT;

      report.passive = false;
      report.violations.clear();

      if (model.rows != model.cols || model.rows < 1)
      {
        return -1;
      }

      RealT low  = RealT{1};
      RealT high = RealT{1};
      bool  scaled = false;
      for (const auto& pole : model.poles)
      {
        const RealT scale = std::abs(pole);
        if (scale <= RealT{0})
        {
          continue;
        }
        low    = scaled ? std::min(low, scale) : scale;
        high   = scaled ? std::max(high, scale) : scale;
        scaled = true;
      }
      low  = low / RealT{SCREEN_RANGE_PAD};
      high = high * RealT{SCREEN_RANGE_PAD};

      std::vector<RealT> grid(SCREEN_POINTS, RealT{0});
      std::vector<RealT> floors(SCREEN_POINTS, RealT{0});
      for (size_t k = 0; k < SCREEN_POINTS; ++k)
      {
        const RealT fraction = static_cast<RealT>(k) /
                               static_cast<RealT>(SCREEN_POINTS - 1);
        grid[k] = std::exp(std::log(low) +
                           fraction * (std::log(high) - std::log(low)));
        if (Detail::conductanceFloor(model, grid[k], floors[k]) != 0)
        {
          return -2;
        }
      }

      // Refine a sign change of the conductance floor by bisection.
      const auto crossing = [&model](RealT inside, RealT outside) -> RealT
      {
        for (size_t step = 0; step < BISECTION_STEPS; ++step)
        {
          const RealT middle = std::sqrt(inside * outside);
          RealT       floor  = RealT{0};
          if (Detail::conductanceFloor(model, middle, floor) != 0)
          {
            break;
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
        return std::sqrt(inside * outside);
      };

      using ReportT   = PassivityReport<scalar_type, index_type>;
      using Violation = typename ReportT::Violation;

      bool      violating = floors[0] < RealT{0};
      Violation band{grid[0], grid[0]};

      for (size_t k = 1; k < SCREEN_POINTS; ++k)
      {
        const bool now = floors[k] < RealT{0};
        if (now && !violating)
        {
          band.omega_start = crossing(grid[k], grid[k - 1]);
          violating        = true;
        }
        else if (!now && violating)
        {
          band.omega_end = crossing(grid[k - 1], grid[k]);
          report.violations.push_back(band);
          violating = false;
        }
      }
      if (violating)
      {
        band.omega_end = grid[SCREEN_POINTS - 1];
        report.violations.push_back(band);
      }

      report.passive = report.violations.empty();
      return 0;
    }

    template int
    assessPassivity<double, long int>(const RationalModel<double, long int>&,
                                      PassivityReport<double, long int>&);
    template int
    assessPassivity<double, size_t>(const RationalModel<double, size_t>&,
                                    PassivityReport<double, size_t>&);
  } // namespace Optimization
} // namespace GridKit
```
