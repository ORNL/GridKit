#pragma once

#include <vector>

#include <GridKit/ScalarTraits.hpp>

namespace AnalysisManager
{
  enum class OptimizationStatus
  {
    NOT_RUN,
    SOLVED,
    ACCEPTABLE,
    INFEASIBLE,
    LIMIT_REACHED,
    USER_STOPPED,
    INVALID_PROBLEM,
    SOLVER_ERROR
  };

  template <class ScalarT, typename IdxT>
  struct OptimizationResult
  {
    using RealT = typename GridKit::ScalarTraits<ScalarT>::RealT;

    OptimizationStatus status{OptimizationStatus::NOT_RUN};
    RealT              objective{};
    RealT              constraint_violation{};
    IdxT               iterations{};

    std::vector<ScalarT> primal;
    std::vector<RealT>   constraints;
    std::vector<RealT>   constraint_duals;
    std::vector<RealT>   lower_bound_duals;
    std::vector<RealT>   upper_bound_duals;

    bool solved() const
    {
      return status == OptimizationStatus::SOLVED || status == OptimizationStatus::ACCEPTABLE;
    }
  };

} // namespace AnalysisManager
