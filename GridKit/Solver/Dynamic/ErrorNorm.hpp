#pragma once

#include <GridKit/LinearAlgebra/Vector/Vector.hpp>
#include <GridKit/LinearAlgebra/Vector/VectorHandler.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace Integrator
{
  /**
   * @brief Interface for error norms. Used to calculate the `err` parameter in `StepController::nextStep` based on a residual state error vector.
   *
   */
  template <class ScalarT, typename IdxT>
  class ErrorNorm
  {
    using State = GridKit::LinearAlgebra::Vector<ScalarT, IdxT>;
    using RealT = typename GridKit::ScalarTraits<ScalarT>::RealT;

  public:
    /**
     * @brief Calculate an error to be used by a step controller. Typically, an error > 1 indicates an error which does not meet tolerances, while
     *        an error < 1 indicates an error which meets tolerances. For that reason, tolerances should be included in the calculation of the error.
     *
     * @param err The state error residual being measured.
     * @param y The state that the error was calculated from. Can be used for proper relative error normalization.
     * @param yprev The state from the previous step. Can be used for proper relative error normalization.
     * @param handler A vector handler which can be used to facilitate vector operations.
     * @param memspace The memory space which vector operations should be performed in/.
     * @return The error.
     *
     * @todo Allow this method to fail, since it will likely involve linear algebra calls.
     */
    virtual RealT errorNorm(State& err, State& y, State& yprev, GridKit::LinearAlgebra::VectorHandler<ScalarT, IdxT>& handler, GridKit::memory::MemorySpace memspace) const = 0;
  };
} // namespace Integrator