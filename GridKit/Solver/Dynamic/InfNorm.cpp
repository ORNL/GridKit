#include "InfNorm.hpp"

#include <format>

#include <GridKit/Constants.hpp>

namespace Integrator
{
  /**
   * @brief Calculate the infinity error norm as
   *
   * \f[e = \max_i\frac{|\hat{e}_i|}{Atol_i + Rtol \cdot \max\{|y_{0i}|, |y_{1i}|\}},\f]
   *
   * where \f(y_0\f) is the initial state, \f(y_1\f) is the next state, and \f(\hat{e}\f) is the estimated error made in calculating
   * the next state (typically \f(\hat{e} = y_1 - \hat{y}_1\f) for some different-order approximation \f(\hat{y}_1\f)).
   *
   * @param err \f(\hat{e}\f) in the above formula.
   * @param y \f(y_1\f) in the above formula.
   * @param yprev \f(y_0\f) in the above formula.
   * @param handler The handler to be used for performing linear algebra operations.
   * @param memspace The memory space to be used for performing linear lagebra operations.
   * @see `Rosenbrock::errorEstimate()`
   */
  template <class ScalarT, typename IdxT>
  InfNorm<ScalarT, IdxT>::RealT InfNorm<ScalarT, IdxT>::errorNorm(State& err, State& y, State& yprev, GridKit::LinearAlgebra::VectorHandler<ScalarT, IdxT>& handler, GridKit::memory::MemorySpace memspace) const
  {
    if (int err_code = workspace_.out_->copyFromExternal(&err, memspace, memspace))
    {
      throw std::format("GridKit::LinearAlgebra::Vector::copyFromExternal failed with error code {}", err_code);
    }
    if (int err_code = workspace_.scale_->copyFromExternal(&y, memspace, memspace))
    {
      throw std::format("GridKit::LinearAlgebra::Vector::copyFromExternal failed with error code {}", err_code);
    }
    if (int err_code = workspace_.yprev_abs_->copyFromExternal(&yprev, memspace, memspace))
    {
      throw std::format("GridKit::LinearAlgebra::Vector::copyFromExternal failed with error code {}", err_code);
    }

    handler.abs(workspace_.scale_.get(), workspace_.scale_.get(), memspace);
    handler.abs(workspace_.yprev_abs_.get(), workspace_.scale_.get(), memspace);
    handler.max(workspace_.yprev_abs_.get(), workspace_.scale_.get(), workspace_.scale_.get(), memspace);

    // TODO: This scal shouldn't be necessary, but axpy doesn't support scaling the y parameter. In the future,
    // the scaling should be able to be put on the next axpy.
    handler.scal(params_.rel_tol_, workspace_.scale_.get(), memspace);
    handler.axpy(GridKit::ONE<RealT>, params_.abs_tol_.get(), workspace_.scale_.get(), memspace);
    handler.diagSolve(workspace_.scale_.get(), workspace_.out_.get(), memspace);

    return handler.amax(workspace_.out_.get(), memspace);
  }
} // namespace Integrator
