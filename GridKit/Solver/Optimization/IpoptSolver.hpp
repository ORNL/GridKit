#pragma once

#include <string>

#include <GridKit/Solver/Optimization/OptimizationResult.hpp>
#include <GridKit/Solver/Optimization/OptimizationSolver.hpp>

namespace AnalysisManager
{
  namespace IpoptInterface
  {
    template <class ScalarT, typename IdxT>
    struct Options
    {
      using RealT = typename GridKit::ScalarTraits<ScalarT>::RealT;

      RealT       tolerance{1.0e-8};
      int         maximum_iterations{500};
      int         print_level{0};
      std::string linear_solver;
    };

    /**
     * @brief Ipopt optimization solver for a GridKit model evaluator.
     *
     * The model remains independent of Ipopt. All TNLP callback translation is
     * private to this solver implementation.
     *
     * @pre The model is allocated and initialized. Its variable, residual, and
     *      derivative storage is sized, and its CSR sparsity pattern is fixed.
     */
    template <class ScalarT, typename IdxT>
    class Solver : public OptimizationSolver<ScalarT, IdxT>
    {
      using OptimizationSolver<ScalarT, IdxT>::model_;

    public:
      using EvaluatorT = GridKit::Model::Evaluator<ScalarT, IdxT>;
      using OptionsT   = Options<ScalarT, IdxT>;
      using ResultT    = OptimizationResult<ScalarT, IdxT>;

      explicit Solver(EvaluatorT* model);
      ~Solver() override = default;

      ResultT solve(const OptionsT& options = {});

    private:
      class Callbacks;
    };

  } // namespace IpoptInterface
} // namespace AnalysisManager
