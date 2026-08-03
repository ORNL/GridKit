```cpp
/**
 * @file IpoptSolver.hpp
 *
 * @brief Interior-point solver backend wrapping Ipopt.
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#pragma once

#include <string>
#include <vector>

#include <GridKit/Solver/Optimization/OptimizationModel.hpp>
#include <GridKit/Solver/Optimization/OptimizationSolver.hpp>

namespace GridKit
{
  namespace Optimization
  {
    /// Hessian supply strategy for the interior-point iteration.
    enum class HessianMode
    {
      EXACT,          ///< The model provides the Lagrangian Hessian.
      LIMITED_MEMORY  ///< Quasi-Newton approximation inside Ipopt.
    };

    /**
     * @brief Solves an OptimizationModel with Ipopt, mirroring how Ida
     *        wraps SUNDIALS for the Dynamic solver family.
     *
     * The public interface exposes no Ipopt types; the adapter lives
     * entirely in the implementation, which builds only under
     * GRIDKIT_ENABLE_IPOPT.
     */
    template <typename scalar_type, typename index_type>
    class IpoptSolver final : public OptimizationSolver<scalar_type,
                                                        index_type>
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename GridKit::ScalarTraits<ScalarT>::RealT;

      struct Parameters
      {
        RealT       tolerance      = 1.0e-8;
        IdxT        max_iterations = 200;
        int         print_level    = 0;
        HessianMode hessian        = HessianMode::EXACT;
      };

      struct Stats
      {
        IdxT  iterations{0};
        RealT final_objective{0.0};
        int   status{0}; ///< Raw Ipopt application return status.

        std::string report() const;
      };

      explicit IpoptSolver(OptimizationModel<ScalarT, IdxT>* model,
                           Parameters params = {});

      /**
       * @brief Solve the bound model.
       *
       * @return 0 at an optimal point, 1 at an acceptable point, 2 when
       *         the iteration limit was reached, negative on a hard
       *         failure
       */
      [[nodiscard("May fail. Check error code.")]]
      int solve(std::vector<RealT>& x) override;

      const Stats& getStats() const
      {
        return stats_;
      }

    private:
      Parameters params_;
      Stats      stats_;
    };
  } // namespace Optimization
} // namespace GridKit
```
