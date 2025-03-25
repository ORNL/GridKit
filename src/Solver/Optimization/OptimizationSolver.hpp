
#ifndef _OPTIMIZATION_SOLVER_HPP_
#define _OPTIMIZATION_SOLVER_HPP_

#include "Model/Evaluator.hpp"
#include <Solver/Dynamic/Ida.hpp>

namespace AnalysisManager
{
  template <class ScalarT, typename IdxT>
  class DynamicSolver;

  template <class ScalarT, typename IdxT>
  class OptimizationSolver
  {
  public:
    OptimizationSolver()
    {
    }

    OptimizationSolver(Sundials::Ida<ScalarT, IdxT>* integrator)
      : integrator_(integrator)
    {
    }

    virtual ~OptimizationSolver()
    {
    }

  protected:
    GridKit::Model::Evaluator<ScalarT, IdxT>* model_;
    Sundials::Ida<ScalarT, IdxT>*             integrator_;
  };

} // namespace AnalysisManager

#endif // _OPTIMIZATION_SOLVER_HPP_
