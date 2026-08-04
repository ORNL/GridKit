#pragma once

#include <IpTNLP.hpp>

#include <GridKit/Optimization/Evaluator.hpp>

namespace GridKit
{
  namespace Optimization
  {
    /**
     * @brief Exact-derivative Ipopt adapter for an Optimization::Evaluator.
     *
     * @pre The evaluator has been allocated and initialized.
     * @pre Exact Jacobian and Hessian capabilities are available.
     * @pre The evaluator is not reallocated and its topology does not change
     *      while the adapter exists.
     */
    template <class scalar_type, typename index_type>
    class IpoptAdapter final : public Ipopt::TNLP
    {
    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using EvaluatorT = Evaluator<ScalarT, IdxT>;
      using RealT      = typename EvaluatorT::RealT;
      using CsrMatrixT = typename EvaluatorT::CsrMatrixT;

      using Index                     = Ipopt::Index;
      using Number                    = Ipopt::Number;
      using SolverReturn              = Ipopt::SolverReturn;
      using IpoptCalculatedQuantities = Ipopt::IpoptCalculatedQuantities;
      using IpoptData                 = Ipopt::IpoptData;

      explicit IpoptAdapter(EvaluatorT& model);
      ~IpoptAdapter() override = default;

      bool get_nlp_info(Index&          variable_count,
                        Index&          constraint_count,
                        Index&          jacobian_nonzeros,
                        Index&          hessian_nonzeros,
                        IndexStyleEnum& index_style) override;

      bool get_bounds_info(Index   variable_count,
                           Number* variable_lower_bounds,
                           Number* variable_upper_bounds,
                           Index   constraint_count,
                           Number* constraint_lower_bounds,
                           Number* constraint_upper_bounds) override;

      bool get_starting_point(Index   variable_count,
                              bool    initialize_variables,
                              Number* variables,
                              bool    initialize_bound_multipliers,
                              Number* lower_bound_multipliers,
                              Number* upper_bound_multipliers,
                              Index   constraint_count,
                              bool    initialize_constraint_multipliers,
                              Number* constraint_multipliers) override;

      bool eval_f(Index         variable_count,
                  const Number* variables,
                  bool          new_variables,
                  Number&       objective) override;

      bool eval_grad_f(Index         variable_count,
                       const Number* variables,
                       bool          new_variables,
                       Number*       gradient) override;

      bool eval_g(Index         variable_count,
                  const Number* variables,
                  bool          new_variables,
                  Index         constraint_count,
                  Number*       constraints) override;

      bool eval_jac_g(Index         variable_count,
                      const Number* variables,
                      bool          new_variables,
                      Index         constraint_count,
                      Index         nonzero_count,
                      Index*        rows,
                      Index*        columns,
                      Number*       values) override;

      bool eval_h(Index         variable_count,
                  const Number* variables,
                  bool          new_variables,
                  Number        objective_factor,
                  Index         constraint_count,
                  const Number* constraint_multipliers,
                  bool          new_constraint_multipliers,
                  Index         nonzero_count,
                  Index*        rows,
                  Index*        columns,
                  Number*       values) override;

      void finalize_solution(SolverReturn               status,
                             Index                      variable_count,
                             const Number*              variables,
                             const Number*              lower_bound_multipliers,
                             const Number*              upper_bound_multipliers,
                             Index                      constraint_count,
                             const Number*              constraints,
                             const Number*              constraint_multipliers,
                             Number                     objective,
                             const IpoptData*           ipopt_data,
                             IpoptCalculatedQuantities* ipopt_quantities) override;

    private:
      bool syncVariables(Index variable_count, const Number* variables);

      EvaluatorT& model_;
      Index       variable_count_{0};
      Index       constraint_count_{0};
      Index       jacobian_nonzeros_{0};
      Index       hessian_nonzeros_{0};
    };
  } // namespace Optimization
} // namespace GridKit
