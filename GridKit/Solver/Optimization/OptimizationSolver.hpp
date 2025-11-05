
#pragma once

#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>

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
