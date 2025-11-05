
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

  protected:
    GridKit::Model::Evaluator<ScalarT, IdxT>* model_;
  };

} // namespace AnalysisManager
