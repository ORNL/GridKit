
#pragma once

#include <IpTNLP.hpp>

#include <GridKit/Solver/Optimization/OptimizationSolver.hpp>

namespace AnalysisManager
{

  namespace IpoptInterface
  {

    /**
     * Implementation of Ipopt's pure virtual TNLP class.
     *
     * TNLP defines Ipopt's interface to the model. This is in fact
     * the model evaluator interface to Ipopt. In this case however,
     * the model evaluator calls dynamic solver to compute the objective
     * and the gradient.
     *
     */
    template <class ScalarT, typename IdxT>
    class DynamicObjective : public Ipopt::TNLP, public OptimizationSolver<ScalarT, IdxT>
    {
      using OptimizationSolver<ScalarT, IdxT>::integrator_;
      using OptimizationSolver<ScalarT, IdxT>::model_;

      typedef typename GridKit::ScalarTraits<ScalarT>::RealT RealT;

      typedef Ipopt::Index                     Index;
      typedef Ipopt::Number                    Number;
      typedef Ipopt::SolverReturn              SolverReturn;
      typedef Ipopt::IpoptCalculatedQuantities IpoptCalculatedQuantities;
      typedef Ipopt::IpoptData                 IpoptData;

    public:
      DynamicObjective(Sundials::Ida<ScalarT, IdxT>* integrator);
      virtual ~DynamicObjective();

      /// Returns sizes of the model components
      virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style);

      /// Returns problem bounds
      virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u);

      /// Initialize optimization
      virtual bool get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda);

      /// Evaluate objective
      virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

      /// Evaluate objective gradient
      virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

      /// Evaluate constraint residuals (not used here)
      virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

      /// Evaluate Jacobian (not used here)
      virtual bool eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nele_jac, Index* iRow, Index* jCol, Number* values);

      /// Evaluate Hessian (have Ipopt estimate Hessian)
      virtual bool eval_h(Index n, const Number* x, bool new_x, Number obj_factor, Index m, const Number* lambda, bool new_lambda, Index nele_hess, Index* iRow, Index* jCol, Number* values);

      /// Postprocessing of the results (not used here)
      virtual void finalize_solution(SolverReturn               status,
                                     Index                      n,
                                     const Number*              x,
                                     const Number*              z_L,
                                     const Number*              z_U,
                                     Index                      m,
                                     const Number*              g,
                                     const Number*              lambda,
                                     Number                     obj_value,
                                     const IpoptData*           ip_data,
                                     IpoptCalculatedQuantities* ip_cq);

    private:
      RealT t_init_;
      RealT t_final_;
      int       nout_;
    };

  } // namespace IpoptInterface
} // namespace AnalysisManager
