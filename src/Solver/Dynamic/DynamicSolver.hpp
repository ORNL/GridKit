
#ifndef _DYNAMIC_SOLVER_HPP_
#define _DYNAMIC_SOLVER_HPP_

#include "Model/Evaluator.hpp"

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

#endif // _DYNAMIC_SOLVER_HPP_
