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
     * @brief Implementation of Ipopt's pure virtual TNLP class for a
     * nonlinear program with exact sparse derivatives
     *
     * Ipopt sees the constraints \f$g = [c; d]\f$, the equality constraints
     * followed by the inequality constraints, with the bounds \f$[0; d_l]\f$
     * and \f$[0; d_u]\f$. The constraint Jacobians and the lower triangle of
     * the Lagrangian Hessian come from the model in CSR format, with the
     * patterns set by the model's `allocate()`.
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
      /**
       * @param[in] model - Allocated and initialized model, which holds the
       * solution after the solve
       */
      explicit IpoptInterface(ModelT* model);

      /// Sizes of the variables, the constraints, and the derivative nonzeros, with C-style indices
      bool get_nlp_info(ipindex& n, ipindex& m, ipindex& nnz_jac_g, ipindex& nnz_h_lag, IndexStyleEnum& index_style) override;

      /// Variable bounds, then zero bounds on the equality constraints followed by the inequality constraint bounds
      bool get_bounds_info(ipindex n, ipnumber* x_l, ipnumber* x_u, ipindex m, ipnumber* g_l, ipnumber* g_u) override;

      /// Variables from the model's `x()`, with no initial multipliers
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

      /// Equality constraints followed by the inequality constraints, \f$g = [c; d]\f$
      bool eval_g(ipindex n, const ipnumber* x, bool new_x, ipindex m, ipnumber* g) override;

      /// Equality constraint Jacobian rows followed by the inequality ones, or their structure when `values` is null
      bool eval_jac_g(ipindex         n,
                      const ipnumber* x,
                      bool            new_x,
                      ipindex         m,
                      ipindex         nele_jac,
                      ipindex*        iRow,
                      ipindex*        jCol,
                      ipnumber*       values) override;

      /// Lower Lagrangian Hessian with the equality multipliers followed by the inequality ones, or its structure when `values` is null
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

      /// Leave the final point and its objective in the model
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
      /// Copy `x` into the model variables
      void setVariables(const ipnumber* x);

      /// Coordinates of the entries of `matrix` in storage order, with the rows shifted by `row_offset`
      static void structure(CsrMatrixT& matrix, ipindex row_offset, ipindex* iRow, ipindex* jCol);

      ModelT* model_;
    };
  } // namespace SteadyState
} // namespace AnalysisManager
