#pragma once

#include <memory>

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
      using State = GridKit::LinearAlgebra::Vector<ScalarT, IdxT>;
      using RealT = ErrorNorm<ScalarT, IdxT>::RealT;

    public:
      /**
       * @brief Configurable error tolerances.
       */
      struct Parameters
      {
        /**
         * @brief Component-wise absolute tolerances.
         */
        std::unique_ptr<State> abs_tol_;
        /**
         * @brief Relative tolerance applied to the larger magnitude of the current and previous states.
         */
        RealT                  rel_tol_;
      } params_;

      /**
       * @brief Construct an RMS norm that owns the supplied absolute-tolerance vector.
       *
       * @param params Error tolerances.
       */
      explicit RmsNorm(Parameters&& params)
        : params_(std::move(params))
      {
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
