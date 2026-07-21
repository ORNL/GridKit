
#pragma once

#include <GridKit/Model/Evaluator.hpp>

namespace AnalysisManager
{
  template <typename scalar_type, typename index_type>
  class DynamicSolver
  {
  public:
    using ScalarT = scalar_type;
    using IdxT    = index_type;

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

    void setTolerance(ScalarT rel_tol)
    {
      setTolerance(rel_tol, 0);
    }

    virtual void setTolerance(ScalarT rel_tol, ScalarT abs_tol_override) = 0;
    virtual void setMaxSteps(IdxT msa)                                   = 0;

  protected:
    GridKit::Model::Evaluator<ScalarT, IdxT>* model_;
  };

} // namespace AnalysisManager
