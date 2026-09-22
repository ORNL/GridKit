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
     */
    template <class ScalarT, typename IdxT>
    InfNorm<ScalarT, IdxT>::RealT InfNorm<ScalarT, IdxT>::errorNorm(State& err, State& y, State& yprev, GridKit::LinearAlgebra::VectorHandler<ScalarT, IdxT>& handler, GridKit::memory::MemorySpace memspace) const
    {
      return handler.weightedInfNorm(&err, &y, &yprev, params_.abs_tol_.get(), params_.rel_tol_, memspace);
    }

    template class InfNorm<double, int>;
  } // namespace NativeDynamicSolver
} // namespace AnalysisManager
