
#include "DynamicConstraint.hpp"

#include <cassert>
#include <iostream>

namespace AnalysisManager
{
  namespace IpoptInterface
  {

    template <class ScalarT, typename IdxT>
    DynamicConstraint<ScalarT, IdxT>::DynamicConstraint(Sundials::Ida<ScalarT, IdxT>* integrator)
      : OptimizationSolver<ScalarT, IdxT>(integrator),
        t_init_(integrator_->getInitialTime()),
        t_final_(integrator_->getFinalTime()),
        nout_(integrator_->getNumberOutputTimes())
    {
      model_ = integrator_->getModel();
    }

    template <class ScalarT, typename IdxT>
    DynamicConstraint<ScalarT, IdxT>::~DynamicConstraint()
    {
    }

    template <class ScalarT, typename IdxT>
    bool DynamicConstraint<ScalarT, IdxT>::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style)
    {
      // This code handles one objective function
      assert(model_->sizeQuadrature() == 1);

      // Number of parameters is size of the system plus 1 fictitious parameter
      // to store the objective value.
      n = static_cast<Index>(model_->sizeParams() + 1);

      // There is one constraint
      m = 1;

      // Jacobian is a dense row matrix of length n+1.
      nnz_jac_g = static_cast<Index>(model_->sizeParams() + 1);

      // Using numerical Hessian.
      nnz_h_lag = 0;

      // Use the C index style (0-based) for row/column entries
      index_style = C_STYLE;

      return true;
    }

    template <class ScalarT, typename IdxT>
    bool DynamicConstraint<ScalarT, IdxT>::get_bounds_info([[maybe_unused]] Index n,
                                                           Number*                x_l,
                                                           Number*                x_u,
                                                           [[maybe_unused]] Index m,
                                                           Number*                g_l,
                                                           Number*                g_u)
    {
      // Check if sizes are set correctly
      assert(n == (Index) (model_->sizeParams() + 1));
      assert(m == 1);

      // Get boundaries for the optimization parameters
      const auto* param_lo = model_->param_lo().getData();
      const auto* param_up = model_->param_up().getData();
      for (IdxT i = 0; i < model_->sizeParams(); ++i)
      {
        x_l[i] = param_lo[static_cast<size_t>(i)];
        x_u[i] = param_up[static_cast<size_t>(i)];
      }

      // No boundaries for fictitious parameter x[n]
      x_l[model_->sizeParams()] = -1e20;
      x_u[model_->sizeParams()] = +1e20;

      // Set constraint g[0] to be equality constraint g[0] = 0
      g_l[0] = 0.0;
      g_u[0] = 0.0;

      return true;
    }

    template <class ScalarT, typename IdxT>
    bool DynamicConstraint<ScalarT, IdxT>::get_starting_point([[maybe_unused]] Index   n,
                                                              [[maybe_unused]] bool    init_x,
                                                              Number*                  x,
                                                              [[maybe_unused]] bool    init_z,
                                                              [[maybe_unused]] Number* z_L,
                                                              [[maybe_unused]] Number* z_U,
                                                              [[maybe_unused]] Index   m,
                                                              [[maybe_unused]] bool    init_lambda,
                                                              [[maybe_unused]] Number* lambda)
    {
      // Only initial values for x provided.
      assert(init_x == true);
      assert(init_z == false);
      assert(init_lambda == false);

      // Initialize optimization parameters x
      const auto* param = model_->param().getData();
      for (IdxT i = 0; i < model_->sizeParams(); ++i)
        x[i] = param[static_cast<size_t>(i)];

      // Initialize fictitious parameter x[n-1] to zero
      x[model_->sizeParams()] = 0.0;

      return true;
    }

    template <class ScalarT, typename IdxT>
    bool DynamicConstraint<ScalarT, IdxT>::eval_f([[maybe_unused]] Index n,
                                                  const Number*          x,
                                                  [[maybe_unused]] bool  new_x,
                                                  Number&                obj_value)
    {
      // Set objective to fictitious optimization parameter x[n-1]
      obj_value = x[model_->sizeParams()];

      return true;
    }

    template <class ScalarT, typename IdxT>
    bool DynamicConstraint<ScalarT, IdxT>::eval_grad_f([[maybe_unused]] Index         n,
                                                       [[maybe_unused]] const Number* x,
                                                       [[maybe_unused]] bool          new_x,
                                                       Number*                        grad_f)
    {
      // Objective function equals to the fictitious parameter x[n-1].
      // Gradient, then assumes the simple form:
      for (IdxT i = 0; i < model_->sizeParams(); ++i)
        grad_f[i] = 0.0;
      grad_f[model_->sizeParams()] = 1.0;

      return true;
    }

    template <class ScalarT, typename IdxT>
    bool DynamicConstraint<ScalarT, IdxT>::eval_g([[maybe_unused]] Index n,
                                                  const Number*          x,
                                                  [[maybe_unused]] bool  new_x,
                                                  [[maybe_unused]] Index m,
                                                  Number*                g)
    {
      // Update optimization parameters
      auto* param = model_->param().getData();
      for (IdxT i = 0; i < model_->sizeParams(); ++i)
      {
        param[static_cast<size_t>(i)] = x[i];
        // std::cout << "x[" << i << "] = " << x[i] << "\n";
      }
      model_->param().setDataUpdated();

      // Evaluate objective function
      integrator_->getSavedInitialCondition();
      integrator_->initializeSimulation(t_init_);
      integrator_->initializeQuadrature();

      int status = 0;
      status     = integrator_->runSimulationQuadrature(t_final_, nout_);
      if (status)
      {
        std::cerr << "Integration failed when using Pm = " << x[0] << "\n";
        return false;
      }

      // For now assumes only one forward integrand and multiple optimization parameters.
      g[0] = (integrator_->getIntegral())[0] - x[model_->sizeParams()];
      // std::cout << "constraint:" << g[0] << std::endl;
      return true;
    }

    template <class ScalarT, typename IdxT>
    bool DynamicConstraint<ScalarT, IdxT>::eval_jac_g([[maybe_unused]] Index n,
                                                      const Number*          x,
                                                      [[maybe_unused]] bool  new_x,
                                                      [[maybe_unused]] Index m,
                                                      [[maybe_unused]] Index nele_jac,
                                                      Index*                 iRow,
                                                      Index*                 jCol,
                                                      Number*                values)
    {
      // Set Jacobian sparsity pattern ...
      if (!values)
      {
        for (IdxT i = 0; i < model_->sizeParams(); ++i)
        {
          iRow[i] = 0;
          jCol[i] = static_cast<Index>(i);
        }
        iRow[model_->sizeParams()] = 0;
        jCol[model_->sizeParams()] = static_cast<Index>(model_->sizeParams());
      }
      // ... or compute Jacobian derivatives
      else
      {
        // Update optimization parameters
        auto* param = model_->param().getData();
        for (IdxT i = 0; i < model_->sizeParams(); ++i)
          param[static_cast<size_t>(i)] = x[i];
        model_->param().setDataUpdated();

        // evaluate the gradient of the objective function grad_{x} f(x)
        // This is creating and deleting adjoint system for each iteration!
        // Currently there is no more efficient solution.
        integrator_->initializeAdjoint();

        integrator_->getSavedInitialCondition();
        integrator_->initializeSimulation(t_init_);
        integrator_->initializeQuadrature();

        int status = 0;
        status     = integrator_->runForwardSimulation(t_final_, nout_);
        if (status)
        {
          std::cerr << "Forward integration for adjoint solution failed when using Pm = " << x[0] << "\n";
          return false;
        }

        integrator_->initializeBackwardSimulation(t_final_);

        status = integrator_->runBackwardSimulation(t_init_);
        if (status)
        {
          std::cerr << "Backward integration for adjoint solution failed when using Pm = " << x[0] << "\n";
          return false;
        }

        // For now assumes only one forward integrand and multiple optimization parameters.
        for (IdxT i = 0; i < model_->sizeParams(); ++i)
        {
          values[i] = -((integrator_->getAdjointIntegral())[i]);
        }
        values[model_->sizeParams()] = -1.0;

        integrator_->deleteAdjoint();
      }
      return true;
    }

    template <class ScalarT, typename IdxT>
    bool DynamicConstraint<ScalarT, IdxT>::eval_h([[maybe_unused]] Index         n,
                                                  [[maybe_unused]] const Number* x,
                                                  [[maybe_unused]] bool          new_x,
                                                  [[maybe_unused]] Number        obj_factor,
                                                  [[maybe_unused]] Index         m,
                                                  [[maybe_unused]] const Number* lambda,
                                                  [[maybe_unused]] bool          new_lambda,
                                                  [[maybe_unused]] Index         nele_hess,
                                                  [[maybe_unused]] Index*        iRow,
                                                  [[maybe_unused]] Index*        jCol,
                                                  [[maybe_unused]] Number*       values)
    {
      return true;
    }

    template <class ScalarT, typename IdxT>
    void DynamicConstraint<ScalarT, IdxT>::finalize_solution([[maybe_unused]] SolverReturn               status,
                                                             [[maybe_unused]] Index                      n,
                                                             [[maybe_unused]] const Number*              x,
                                                             [[maybe_unused]] const Number*              z_L,
                                                             [[maybe_unused]] const Number*              z_U,
                                                             [[maybe_unused]] Index                      m,
                                                             [[maybe_unused]] const Number*              g,
                                                             [[maybe_unused]] const Number*              lambda,
                                                             [[maybe_unused]] Number                     obj_value,
                                                             [[maybe_unused]] const IpoptData*           ip_data,
                                                             [[maybe_unused]] IpoptCalculatedQuantities* ip_cq)
    {
    }

    template class DynamicConstraint<sunrealtype, long int>;
    template class DynamicConstraint<sunrealtype, size_t>;

  } // namespace IpoptInterface
} // namespace AnalysisManager
