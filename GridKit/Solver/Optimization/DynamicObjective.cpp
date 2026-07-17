
#include "DynamicObjective.hpp"

#include <cassert>
#include <iostream>

namespace AnalysisManager
{
  namespace IpoptInterface
  {

    template <class ScalarT, typename IdxT>
    DynamicObjective<ScalarT, IdxT>::DynamicObjective(Sundials::Ida<ScalarT, IdxT>* integrator,
                                                      RealT                         t_init,
                                                      RealT                         t_final,
                                                      RealT                         dt_monitor)
      : OptimizationSolver<ScalarT, IdxT>(integrator),
        t_init_(t_init),
        t_final_(t_final),
        dt_monitor_(dt_monitor)
    {
      model_ = integrator_->getModel();
    }

    template <class ScalarT, typename IdxT>
    DynamicObjective<ScalarT, IdxT>::~DynamicObjective()
    {
    }

    template <class ScalarT, typename IdxT>
    bool DynamicObjective<ScalarT, IdxT>::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style)
    {
      // This code handles one objective function
      assert(model_->sizeQuadrature() == 1);

      // Number of optimization variables.
      n = static_cast<Index>(model_->sizeParams());

      // There are no constraints
      m = 0;

      // No constraints, empty Jacobian. This is only temporary.
      nnz_jac_g = 0;

      // Using numerical Hessian.
      nnz_h_lag = 0;

      // Use the C index style (0-based) for row/column entries
      index_style = C_STYLE;

      return true;
    }

    template <class ScalarT, typename IdxT>
    bool DynamicObjective<ScalarT, IdxT>::get_bounds_info([[maybe_unused]] Index   n,
                                                          Number*                  x_l,
                                                          Number*                  x_u,
                                                          [[maybe_unused]] Index   m,
                                                          [[maybe_unused]] Number* g_l,
                                                          [[maybe_unused]] Number* g_u)
    {
      // Check if sizes are set correctly
      assert(n == (Index) model_->sizeParams());
      assert(m == 0);

      // Get boundaries for the optimization parameters
      const auto* param_lo = model_->param_lo().getData();
      const auto* param_up = model_->param_up().getData();
      for (IdxT i = 0; i < model_->sizeParams(); ++i)
      {
        x_l[i] = param_lo[static_cast<size_t>(i)];
        x_u[i] = param_up[static_cast<size_t>(i)];
      }

      return true;
    }

    template <class ScalarT, typename IdxT>
    bool DynamicObjective<ScalarT, IdxT>::get_starting_point([[maybe_unused]] Index   n,
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

      return true;
    }

    template <class ScalarT, typename IdxT>
    bool DynamicObjective<ScalarT, IdxT>::eval_f([[maybe_unused]] Index n,
                                                 const Number*          x,
                                                 [[maybe_unused]] bool  new_x,
                                                 Number&                obj_value)
    {
      // Update optimization parameters
      auto* param = model_->param().getData();
      for (IdxT i = 0; i < model_->sizeParams(); ++i)
        param[static_cast<size_t>(i)] = x[i];
      model_->param().setDataUpdated();

      // Evaluate objective function
      integrator_->getSavedInitialCondition();
      integrator_->initializeSimulation(t_init_);
      integrator_->initializeQuadrature();
      integrator_->runSimulationQuadrature(t_final_, dt_monitor_);

      // Assuming objective function is given as the integral (quadrature) 0
      obj_value = (integrator_->getIntegral())[0];

      return true;
    }

    template <class ScalarT, typename IdxT>
    bool DynamicObjective<ScalarT, IdxT>::eval_grad_f([[maybe_unused]] Index n,
                                                      const Number*          x,
                                                      [[maybe_unused]] bool  new_x,
                                                      Number*                grad_f)
    {
      assert(model_->sizeParams() == static_cast<IdxT>(n));
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
      integrator_->runForwardSimulation(t_final_, dt_monitor_);

      integrator_->initializeBackwardSimulation(t_final_);
      integrator_->runBackwardSimulation(t_init_);

      // For now assumes only one forward integrand and multiple optimization parameters.
      for (IdxT i = 0; i < model_->sizeParams(); ++i)
        grad_f[i] = -((integrator_->getAdjointIntegral())[i]);

      integrator_->deleteAdjoint();

      return true;
    }

    template <class ScalarT, typename IdxT>
    bool DynamicObjective<ScalarT, IdxT>::eval_g([[maybe_unused]] Index         n,
                                                 [[maybe_unused]] const Number* x,
                                                 [[maybe_unused]] bool          new_x,
                                                 [[maybe_unused]] Index         m,
                                                 [[maybe_unused]] Number*       g)
    {
      return true;
    }

    template <class ScalarT, typename IdxT>
    bool DynamicObjective<ScalarT, IdxT>::eval_jac_g([[maybe_unused]] Index         n,
                                                     [[maybe_unused]] const Number* x,
                                                     [[maybe_unused]] bool          new_x,
                                                     [[maybe_unused]] Index         m,
                                                     [[maybe_unused]] Index         nele_jac,
                                                     [[maybe_unused]] Index*        iRow,
                                                     [[maybe_unused]] Index*        jCol,
                                                     [[maybe_unused]] Number*       values)
    {
      return true;
    }

    template <class ScalarT, typename IdxT>
    bool DynamicObjective<ScalarT, IdxT>::eval_h([[maybe_unused]] Index         n,
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
    void DynamicObjective<ScalarT, IdxT>::finalize_solution([[maybe_unused]] SolverReturn               status,
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

    template class DynamicObjective<Ipopt::Number, Ipopt::Index>;
    template class DynamicObjective<Ipopt::Number, std::size_t>;

  } // namespace IpoptInterface
} // namespace AnalysisManager
