#pragma once

#include <GridKit/Model/Evaluator.hpp>

namespace AnalysisManager
{
  template <class ScalarT, typename IdxT>
  class OptimizationSolver
  {
  public:
    using EvaluatorT = GridKit::Model::Evaluator<ScalarT, IdxT>;

    explicit OptimizationSolver(EvaluatorT* model)
      : model_(model)
    {
    }

    virtual ~OptimizationSolver() = default;

    EvaluatorT* getModel()
    {
      return model_;
    }

    const EvaluatorT* getModel() const
    {
      return model_;
    }

  protected:
    EvaluatorT* model_{nullptr};
  };

} // namespace AnalysisManager
