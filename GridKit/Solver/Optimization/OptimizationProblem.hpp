/**
 * @file OptimizationProblem.hpp
 * @brief Ipopt problem over a `Model::OptimizationEvaluator`.
 */

#pragma once

#include <IpTNLP.hpp>

#include <GridKit/Model/OptimizationEvaluator.hpp>

namespace AnalysisManager
{
  namespace IpoptInterface
  {
    /**
     * @brief Implementation of Ipopt's pure virtual TNLP class for a
     * nonlinear program with exact sparse derivatives
     *
     * The constraint Jacobian and the lower triangle of the Lagrangian
     * Hessian come from the model in CSR format, with the pattern set by the
     * model's `allocate()`.
     */
    template <typename scalar_type, typename index_type>
    class OptimizationProblem : public Ipopt::TNLP
    {
    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using ModelT     = GridKit::Model::OptimizationEvaluator<ScalarT, IdxT>;
      using RealT      = typename ModelT::RealT;
      using CsrMatrixT = typename ModelT::CsrMatrixT;

      using Index                     = Ipopt::Index;
      using Number                    = Ipopt::Number;
      using SolverReturn              = Ipopt::SolverReturn;
      using IpoptCalculatedQuantities = Ipopt::IpoptCalculatedQuantities;
      using IpoptData                 = Ipopt::IpoptData;

      /**
       * @param[in] model - Allocated model, which holds the solution after
       * the solve
       */
      explicit OptimizationProblem(ModelT* model);

      bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style) override;

      bool get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u) override;

      bool get_starting_point(Index   n,
                              bool    init_x,
                              Number* x,
                              bool    init_z,
                              Number* z_L,
                              Number* z_U,
                              Index   m,
                              bool    init_lambda,
                              Number* lambda) override;

      bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value) override;

      bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f) override;

      bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) override;

      bool eval_jac_g(Index         n,
                      const Number* x,
                      bool          new_x,
                      Index         m,
                      Index         nele_jac,
                      Index*        iRow,
                      Index*        jCol,
                      Number*       values) override;

      bool eval_h(Index         n,
                  const Number* x,
                  bool          new_x,
                  Number        obj_factor,
                  Index         m,
                  const Number* lambda,
                  bool          new_lambda,
                  Index         nele_hess,
                  Index*        iRow,
                  Index*        jCol,
                  Number*       values) override;

      void finalize_solution(SolverReturn               status,
                             Index                      n,
                             const Number*              x,
                             const Number*              z_L,
                             const Number*              z_U,
                             Index                      m,
                             const Number*              g,
                             const Number*              lambda,
                             Number                     obj_value,
                             const IpoptData*           ip_data,
                             IpoptCalculatedQuantities* ip_cq) override;

    private:
      void setVariables(const Number* x);

      static void structure(CsrMatrixT& matrix, Index* iRow, Index* jCol);

      ModelT* model_;
    };
  } // namespace IpoptInterface
} // namespace AnalysisManager
