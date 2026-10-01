#pragma once

#include <GridKit/Model/Evaluator.hpp>

namespace AnalysisManager
{
  template <class ScalarT, typename IdxT>
  class SteadyStateSolver
  {
  public:
    SteadyStateSolver(GridKit::Model::Evaluator<ScalarT, IdxT>* model)
      : model_(model)
    {
    }

    virtual ~SteadyStateSolver()
    {
    }

    GridKit::Model::Evaluator<ScalarT, IdxT>* getModel()
    {
      return model_;
    }

    virtual void setTolerance(ScalarT tol) = 0;

  protected:
    GridKit::Model::Evaluator<ScalarT, IdxT>* model_;
  };

} // namespace AnalysisManager
