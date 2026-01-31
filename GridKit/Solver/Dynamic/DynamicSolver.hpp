
#pragma once

#include <GridKit/Model/Evaluator.hpp>

namespace AnalysisManager
{
  template <class ScalarT, typename IdxT>
  class DynamicSolver
  {
  public:
    DynamicSolver(GridKit::Model::Evaluator<ScalarT, IdxT>* model)
      : model_(model)
    {
    }

    virtual ~DynamicSolver()
    {
    }

    GridKit::Model::Evaluator<ScalarT, IdxT>* getModel()
    {
      return model_;
    }

    virtual void setTimeStep(ScalarT timeStep) = 0;
    virtual void setTolerance(ScalarT relTol, ScalarT absTolFac=1) = 0;
    virtual void setMaxSteps(IdxT msa) = 0;

  protected:
    GridKit::Model::Evaluator<ScalarT, IdxT>* model_;
  };

} // namespace AnalysisManager
