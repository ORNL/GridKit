/**
 * @file IpoptInterface.hpp
 * @brief Ipopt problem over a `Model::OptimizationEvaluator`.
 */

#pragma once

#include <IpTNLP.hpp>

#include <GridKit/Model/OptimizationEvaluator.hpp>

namespace AnalysisManager
{
  namespace SteadyState
  {
    /**
     * @brief Ipopt problem over a `Model::OptimizationEvaluator`, which holds
     * the solution and its multipliers after the solve
     *
     * Unlike `IpoptInterface::DynamicObjective` and
     * `IpoptInterface::DynamicConstraint`, which integrate a dynamic model to
     * evaluate the problem, it takes the evaluator's functions and exact
     * sparse derivatives.
     *
     * @tparam scalar_type - Scalar of the evaluator, which must be `ipnumber`
     * @tparam index_type - Index of the evaluator
     */
    template <typename scalar_type, typename index_type>
    class IpoptInterface : public Ipopt::TNLP
    {
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using ModelT     = GridKit::Model::OptimizationEvaluator<ScalarT, IdxT>;
      using RealT      = typename ModelT::RealT;
      using CsrMatrixT = typename ModelT::CsrMatrixT;

      using SolverReturn              = Ipopt::SolverReturn;
      using IpoptCalculatedQuantities = Ipopt::IpoptCalculatedQuantities;
      using IpoptData                 = Ipopt::IpoptData;

    public:
      /// Solve `model`, which must be allocated and initialized
      explicit IpoptInterface(ModelT* model);

      /// Sizes and nonzero counts, with C-style indices
      bool get_nlp_info(ipindex& n, ipindex& m, ipindex& nnz_jac_g, ipindex& nnz_h_lag, IndexStyleEnum& index_style) override;

      /// Variable and constraint bounds
      bool get_bounds_info(ipindex n, ipnumber* x_l, ipnumber* x_u, ipindex m, ipnumber* g_l, ipnumber* g_u) override;

      /// Variables from `x()`, without initial multipliers
      bool get_starting_point(ipindex   n,
                              bool      init_x,
                              ipnumber* x,
                              bool      init_z,
                              ipnumber* z_L,
                              ipnumber* z_U,
                              ipindex   m,
                              bool      init_lambda,
                              ipnumber* lambda) override;

      /// Objective \f$f(x)\f$
      bool eval_f(ipindex n, const ipnumber* x, bool new_x, ipnumber& obj_value) override;

      /// Objective gradient \f$\nabla f(x)\f$
      bool eval_grad_f(ipindex n, const ipnumber* x, bool new_x, ipnumber* grad_f) override;

      /// Constraints \f$g(x)\f$
      bool eval_g(ipindex n, const ipnumber* x, bool new_x, ipindex m, ipnumber* g) override;

      /// Constraint Jacobian, or its structure when `values` is null
      bool eval_jac_g(ipindex         n,
                      const ipnumber* x,
                      bool            new_x,
                      ipindex         m,
                      ipindex         nele_jac,
                      ipindex*        iRow,
                      ipindex*        jCol,
                      ipnumber*       values) override;

      /// Lower Lagrangian Hessian, or its structure when `values` is null
      bool eval_h(ipindex         n,
                  const ipnumber* x,
                  bool            new_x,
                  ipnumber        obj_factor,
                  ipindex         m,
                  const ipnumber* lambda,
                  bool            new_lambda,
                  ipindex         nele_hess,
                  ipindex*        iRow,
                  ipindex*        jCol,
                  ipnumber*       values) override;

      /// Leave the final point, its multipliers, and its objective in the model
      void finalize_solution(SolverReturn               status,
                             ipindex                    n,
                             const ipnumber*            x,
                             const ipnumber*            z_L,
                             const ipnumber*            z_U,
                             ipindex                    m,
                             const ipnumber*            g,
                             const ipnumber*            lambda,
                             ipnumber                   obj_value,
                             const IpoptData*           ip_data,
                             IpoptCalculatedQuantities* ip_cq) override;

    private:
      /// Copy `values` from Ipopt into a model vector
      template <typename VectorT>
      static void copyInto(VectorT& vector, const ipnumber* values);

      /// Coordinates of the entries of `matrix` in storage order
      static void structure(CsrMatrixT& matrix, ipindex* iRow, ipindex* jCol);

      ModelT* model_;
    };
  } // namespace SteadyState
} // namespace AnalysisManager
