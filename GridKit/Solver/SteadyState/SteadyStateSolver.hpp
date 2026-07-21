#pragma once

#include <GridKit/Model/Evaluator.hpp>

namespace AnalysisManager
{
  template <typename scalar_type, typename index_type>
  class SteadyStateSolver
  {
  public:
    using ScalarT = scalar_type;
    using IdxT    = index_type;

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
