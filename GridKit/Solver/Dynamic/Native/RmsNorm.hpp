#pragma once

#include <memory>
#include <stdexcept>
#include <utility>

#include <GridKit/LinearAlgebra/Vector/Vector.hpp>
#include <GridKit/Solver/Dynamic/Native/ErrorNorm.hpp>

namespace AnalysisManager
{
  namespace NativeDynamicSolver
  {
    /**
     * @brief Root-mean-square norm of component-wise tolerance-scaled errors.
     */
    template <class ScalarT, typename IdxT>
    class RmsNorm : public ErrorNorm<ScalarT, IdxT>
    {
      using State = typename ErrorNorm<ScalarT, IdxT>::State;
      using RealT = typename ErrorNorm<ScalarT, IdxT>::RealT;

      std::unique_ptr<State> abs_tol_;
      RealT                  rel_tol_;

    public:
      /**
       * @brief Construct an RMS norm that owns the supplied absolute-tolerance vector.
       *
       * @param abs_tol Component-wise absolute tolerances.
       * @param rel_tol Relative tolerance applied to the larger magnitude of the current and previous states.
       */
      RmsNorm(std::unique_ptr<State> abs_tol, RealT rel_tol)
        : abs_tol_(std::move(abs_tol)),
          rel_tol_(rel_tol)
      {
        if (!abs_tol_)
          throw std::invalid_argument("RmsNorm requires an absolute-tolerance vector");
      }

      /**
       * @brief Compute the tolerance-scaled RMS norm through the supplied vector handler.
       */
      RealT errorNorm(State&                                                err,
                      State&                                                y,
                      State&                                                yprev,
                      GridKit::LinearAlgebra::VectorHandler<ScalarT, IdxT>& handler,
                      GridKit::memory::MemorySpace                          memspace) const final;
    };
  } // namespace NativeDynamicSolver
} // namespace AnalysisManager
