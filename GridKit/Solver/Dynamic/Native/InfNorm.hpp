#pragma once

#include <memory>

#include <GridKit/LinearAlgebra/Vector/Vector.hpp>
#include <GridKit/Solver/Dynamic/Native/ErrorNorm.hpp>

namespace Integrator
{

  /**
   * @brief The infinity error norm, which requires the error in every component in the system
   *        to meet tolerance.
   *
   */
  template <class ScalarT, typename IdxT>
  class InfNorm : public ErrorNorm<ScalarT, IdxT>
  {
    using State = GridKit::LinearAlgebra::Vector<ScalarT, IdxT>;
    using RealT = ErrorNorm<ScalarT, IdxT>::RealT;

    /**
     * @brief A workspace for the linear algebra operations required to calculate the norm.
     *
     */
    mutable struct
    {
      /**
       * @brief The final vector which will have its norm taken.
       *
       */
      std::unique_ptr<State> out_;
      /**
       * @brief The vector which will be used to scale the error in accordance with the tolerances.
       *
       */
      std::unique_ptr<State> scale_;
      /**
       * @brief The absolute value of yprev. Used to calculate \ref scale_
       *
       */
      std::unique_ptr<State> yprev_abs_;
    } workspace_;

  public:
    /**
     * @brief The configurable parameters of the error norm.
     *
     */
    struct Parameters
    {
      /**
       * @brief A vector of absolute tolerances for each component. The norm will attempt to reject errors
       *        in each component above this tolerance.
       *
       */
      std::unique_ptr<State> abs_tol_;

      /**
       * @brief The relative tolerance. The norm will attempt to reject any error larger in percentage of
       *        the solution's maximum element than this.
       *
       */
      RealT rel_tol_;
    } params_;

    InfNorm(Parameters&& params)
      : params_(std::move(params))
    {
    }

    RealT errorNorm(State& err, State& y, State& yprev, GridKit::LinearAlgebra::VectorHandler<ScalarT, IdxT>& handler, GridKit::memory::MemorySpace memspace) const final;
  };

} // namespace Integrator
