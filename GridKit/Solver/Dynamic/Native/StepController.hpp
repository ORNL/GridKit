#pragma once

#include <cstdint>

#include <GridKit/Solver/Dynamic/Native/StepControl.hpp>

namespace AnalysisManager
{
  namespace NativeDynamicSolver
  {
    /**
     * @brief Interface for step size controllers. Used by `Rosenbrock` integrators to decide when to accept/reject steps and
     *        what size each step should be.
     *
     * @todo It may be best to have \ref usesError() return a reference to the \ref ErrorNorm that should be used.
     */
    template <typename RealT>
    class StepController
    {
    public:
      /**
       * @brief Decide the control flow for the next step, based on information gathered by the integrator about the current step.
       *
       * @param err The estimated error made by the current step. Only calculated if `usesError()` returns true, otherwise it is assumed this
       *            value is not used.
       * @param prev_step The control flow from before the current step was taken. Can be used to accurately update the step size.
       * @param method_order The order of the method being used.
       * @return StepControl
       */
      virtual StepControl<RealT> nextStep(RealT err, StepControl<RealT> prev_step, uint8_t method_order) = 0;
      /**
       * @brief Return whether or not the `nextStep` method implementation uses the `err` parameter. If `false`, this parameter is not calculated.
       *
       */
      virtual bool               usesError() const                                                       = 0;
    };
  } // namespace NativeDynamicSolver
} // namespace AnalysisManager
