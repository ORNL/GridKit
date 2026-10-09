#include "InfNorm.hpp"

#include <stdexcept>

namespace AnalysisManager
{
  namespace NativeDynamicSolver
  {
    /**
     * @brief Calculate the infinity error norm as
     *
     * \f[e = \max_i\frac{|\hat{e}_i|}{Atol_i + Rtol \cdot \max\{|y_{0i}|, |y_{1i}|\}}.\f]
     *
     * where \f(y_0\f) is the initial state, \f(y_1\f) is the next state, and \f(\hat{e}\f) is the estimated error made in calculating
     * the next state (typically \f(\hat{e} = y_1 - \hat{y}_1\f) for some different-order approximation \f(\hat{y}_1\f)).
     *
     * @param err \f(\hat{e}\f) in the above formula.
     * @param y \f(y_1\f) in the above formula.
     * @param yprev \f(y_0\f) in the above formula.
     * @param handler The handler to be used for performing linear algebra operations.
     * @param memspace The memory space to be used for performing linear algebra operations.
     * @see `Rosenbrock::errorEstimate()`
     */
    template <class ScalarT, typename IdxT>
    typename InfNorm<ScalarT, IdxT>::RealT InfNorm<ScalarT, IdxT>::errorNorm(State& err, State& y, State& yprev, GridKit::LinearAlgebra::VectorHandler<ScalarT, IdxT>& handler, GridKit::memory::MemorySpace memspace) const
    {
      return handler.weightedInfNorm(&err, &y, &yprev, abs_tol_.get(), rel_tol_, memspace);
    }

    template class InfNorm<double, int>;
  } // namespace NativeDynamicSolver
} // namespace AnalysisManager
