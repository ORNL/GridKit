/*
 *
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Slaven Peles <peles2@llnl.gov>.
 * LLNL-CODE-718378.
 * All rights reserved.
 *
 * This file is part of GridKit™. For details, see github.com/LLNL/GridKit
 * Please also read the LICENSE file.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * - Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the disclaimer below.
 * - Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the disclaimer (as noted below) in the
 *   documentation and/or other materials provided with the distribution.
 * - Neither the name of the LLNS/LLNL nor the names of its contributors may
 *   be used to endorse or promote products derived from this software without
 *   specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL
 * SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 * OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISINGIN ANY
 * WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Lawrence Livermore National Laboratory is operated by Lawrence Livermore
 * National Security, LLC, for the U.S. Department of Energy, National
 * Nuclear Security Administration under Contract DE-AC52-07NA27344.
 *
 * This document was prepared as an account of work sponsored by an agency
 * of the United States government. Neither the United States government nor
 * Lawrence Livermore National Security, LLC, nor any of their employees
 * makes any warranty, expressed or implied, or assumes any legal liability
 * or responsibility for the accuracy, completeness, or usefulness of any
 * information, apparatus, product, or process disclosed, or represents that
 * its use would not infringe privately owned rights. Reference herein to
 * any specific commercial product, process, or service by trade name,
 * trademark, manufacturer, or otherwise does not necessarily constitute or
 * imply its endorsement, recommendation, or favoring by the United States
 * government or Lawrence Livermore National Security, LLC. The views and
 * opinions of authors expressed herein do not necessarily state or reflect
 * those of the United States government or Lawrence Livermore National
 * Security, LLC, and shall not be used for advertising or product
 * endorsement purposes.
 *
 */


#include <iostream>
#include <cassert>

#include "DynamicConstraint.hpp"

namespace AnalysisManager {
namespace IpoptInterface {

template <class ScalarT, typename IdxT>
DynamicConstraint<ScalarT, IdxT>::DynamicConstraint(Sundials::Ida<ScalarT, IdxT>* integrator)
  : OptimizationSolver<ScalarT, IdxT>(integrator),
    t_init_(integrator_->getInitialTime()),
    t_final_(integrator_->getFinalTime()),
    nout_(integrator_->getNumberOutputTimes())
{
    model_ =  integrator_->getModel();
}

template <class ScalarT, typename IdxT>
DynamicConstraint<ScalarT, IdxT>::~DynamicConstraint()
{
}


template <class ScalarT, typename IdxT>
bool DynamicConstraint<ScalarT, IdxT>::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                                              Index& nnz_h_lag, IndexStyleEnum& index_style)
{
    // This code handles one objective function
    assert(model_->size_quad() == 1);

    // Number of parameters is size of the system plus 1 fictitious parameter
    // to store the objective value.
    n = model_->size_opt() + 1;

    // There is one constraint
    m = 1;

    // Jacobian is a dense row matrix of length n+1.
    nnz_jac_g = model_->size_opt() + 1;

    // Using numerical Hessian.
    nnz_h_lag = 0;

    // Use the C index style (0-based) for row/column entries
    index_style = C_STYLE;

    return true;
}


template <class ScalarT, typename IdxT>
bool DynamicConstraint<ScalarT, IdxT>::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                                 Index m, Number* g_l, Number* g_u)
{
    // Check if sizes are set correctly
    assert(n == (Index) (model_->size_opt() + 1));
    assert(m == 1);

    // Get boundaries for the optimization parameters
    for(IdxT i = 0; i < model_->size_opt(); ++i)
    {
        x_l[i] = model_->param_lo()[i];
        x_u[i] = model_->param_up()[i];
    }

    // No boundaries for fictitious parameter x[n]
    x_l[model_->size_opt()] = -1e20;
    x_u[model_->size_opt()] = +1e20;

    // Set constraint g[0] to be equality constraint g[0] = 0
    g_l[0] = 0.0;
    g_u[0] = 0.0;

    return true;
}


template <class ScalarT, typename IdxT>
bool DynamicConstraint<ScalarT, IdxT>::get_starting_point(Index n, bool init_x, Number* x,
                                                          bool init_z, Number* z_L, Number* z_U,
                                                          Index m, bool init_lambda,
                                                          Number* lambda)
{
    // Only initial values for x provided.
    assert(init_x == true);
    assert(init_z == false);
    assert(init_lambda == false);

    // Initialize optimization parameters x
    for(IdxT i = 0; i < model_->size_opt(); ++i)
        x[i] = model_->param()[i];

    // Initialize fictitious parameter x[n-1] to zero
    x[model_->size_opt()] = 0.0;

    return true;
}


template <class ScalarT, typename IdxT>
bool DynamicConstraint<ScalarT, IdxT>::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
    // Set objective to fictitious optimization parameter x[n-1]
    obj_value = x[model_->size_opt()];

    return true;
}


template <class ScalarT, typename IdxT>
bool DynamicConstraint<ScalarT, IdxT>::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
    // Objective function equals to the fictitious parameter x[n-1].
    // Gradient, then assumes the simple form:
    for(IdxT i = 0; i < model_->size_opt(); ++i)
        grad_f[i] = 0.0;
    grad_f[model_->size_opt()] = 1.0;

    return true;
}


template <class ScalarT, typename IdxT>
bool DynamicConstraint<ScalarT, IdxT>::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
    // Update optimization parameters
    for(IdxT i = 0; i < model_->size_opt(); ++i)
    {
        model_->param()[i] = x[i];
        //std::cout << "x[" << i << "] = " << x[i] << "\n";
    }

    // Evaluate objective function
    integrator_->getSavedInitialCondition();
    integrator_->initializeSimulation(t_init_);
    integrator_->initializeQuadrature();

    int status = 0;
    status = integrator_->runSimulationQuadrature(t_final_, nout_);
    if (status)
    {
        std::cerr << "Integration failed when using Pm = " << x[0] << "\n";
        return false;
    }

    // For now assumes only one forward integrand and multiple optimization parameters.
    g[0] = (integrator_->getIntegral())[0] - x[model_->size_opt()];
    //std::cout << "constraint:" << g[0] << std::endl;
    return true;
}


template <class ScalarT, typename IdxT>
bool DynamicConstraint<ScalarT, IdxT>::eval_jac_g(Index n, const Number* x, bool new_x,
                                                  Index m, Index nele_jac, Index* iRow, Index *jCol,
                                                  Number* values)
{
    // Set Jacobian sparsity pattern ...
    if(!values)
    {
        for(IdxT i = 0; i < model_->size_opt(); ++i)
        {
            iRow[i] = 0;
            jCol[i] = i;
        }
        iRow[model_->size_opt()] = 0;
        jCol[model_->size_opt()] = model_->size_opt();
    }
    // ... or compute Jacobian derivatives
    else
    {
        // Update optimization parameters
        for(IdxT i = 0; i < model_->size_opt(); ++i)
            model_->param()[i] = x[i];

        // evaluate the gradient of the objective function grad_{x} f(x)
        // This is creating and deleting adjoint system for each iteration!
        // Currently there is no more efficient solution.
        integrator_->initializeAdjoint();

        integrator_->getSavedInitialCondition();
        integrator_->initializeSimulation(t_init_);
        integrator_->initializeQuadrature();

        int status = 0;
        status = integrator_->runForwardSimulation(t_final_, nout_);
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
        for(IdxT i = 0; i < model_->size_opt(); ++i)
        {
            values[i] = -((integrator_->getAdjointIntegral())[i]);
        }
        values[model_->size_opt()] = -1.0;

        integrator_->deleteAdjoint();
    }
    return true;
}


template <class ScalarT, typename IdxT>
bool DynamicConstraint<ScalarT, IdxT>::eval_h(Index n, const Number* x, bool new_x,
                                              Number obj_factor, Index m, const Number* lambda,
                                              bool new_lambda, Index nele_hess, Index* iRow,
                                              Index* jCol, Number* values)
{
    return true;
}


template <class ScalarT, typename IdxT>
void DynamicConstraint<ScalarT, IdxT>::finalize_solution(SolverReturn status,
                                                         Index n, const Number* x, const Number* z_L, const Number* z_U,
                                                         Index m, const Number* g, const Number* lambda,
                                                         Number obj_value,
                                                         const IpoptData* ip_data,
                                                         IpoptCalculatedQuantities* ip_cq)
{

}




template class DynamicConstraint<realtype, long int>;
template class DynamicConstraint<realtype, size_t>;

} // namespace IpoptInterface
} // namespace AnalysisManager
