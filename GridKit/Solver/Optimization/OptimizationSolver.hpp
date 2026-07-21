
#pragma once

#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>

namespace AnalysisManager
{
  template <typename scalar_type, typename index_type>
  class DynamicSolver;

  template <typename scalar_type, typename index_type>
  class OptimizationSolver
  {
  public:
    using ScalarT = scalar_type;
    using IdxT    = index_type;

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
