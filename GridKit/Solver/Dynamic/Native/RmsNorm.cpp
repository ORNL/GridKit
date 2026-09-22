#include "RmsNorm.hpp"

#include <stdexcept>

namespace AnalysisManager
{
  namespace NativeDynamicSolver
  {
    /**
     * @brief Calculate the weighted root-mean-square error norm
     *
     * \f[e = \sqrt{\frac{1}{N}\sum_i\left(\frac{\hat{e}_i}
     * {Atol_i + Rtol \cdot \max\{|y_{0i}|, |y_{1i}|\}}\right)^2}.\f]
     *
     * @param err Estimated error.
     * @param y Current state.
     * @param yprev Previous state.
     * @param handler Vector handler used for the fused reduction.
     * @param memspace Memory space in which to perform the reduction.
     * @return The weighted RMS norm.
     */
    template <class ScalarT, typename IdxT>
    RmsNorm<ScalarT, IdxT>::RealT RmsNorm<ScalarT, IdxT>::errorNorm(State&                                                err,
                                                                    State&                                                y,
                                                                    State&                                                yprev,
                                                                    GridKit::LinearAlgebra::VectorHandler<ScalarT, IdxT>& handler,
                                                                    GridKit::memory::MemorySpace                          memspace) const
    {
      return handler.weightedRmsNorm(&err, &y, &yprev, params_.abs_tol_.get(), params_.rel_tol_, memspace);
    }

    template class RmsNorm<double, int>;
  } // namespace NativeDynamicSolver
} // namespace AnalysisManager
