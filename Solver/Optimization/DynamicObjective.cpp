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

#include "DynamicObjective.hpp"

namespace AnalysisManager {
namespace IpoptInterface {

template <class ScalarT, typename IdxT>
DynamicObjective<ScalarT, IdxT>::DynamicObjective(Sundials::Ida<ScalarT, IdxT>* integrator)
  : OptimizationSolver<ScalarT, IdxT>(integrator),
    t_init_(integrator_->getInitialTime()),
    t_final_(integrator_->getFinalTime()),
    nout_(integrator_->getNumberOutputTimes())
{
    model_ =  integrator_->getModel();
}

template <class ScalarT, typename IdxT>
DynamicObjective<ScalarT, IdxT>::~DynamicObjective()
{
}


template <class ScalarT, typename IdxT>
bool DynamicObjective<ScalarT, IdxT>::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                                                   Index& nnz_h_lag, IndexStyleEnum& index_style)
{
    // This code handles one objective function
    assert(model_->size_quad() == 1);
    
    // Number of optimization variables.
    n = model_->size_opt();

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
bool DynamicObjective<ScalarT, IdxT>::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                                      Index m, Number* g_l, Number* g_u)
{
    // Check if sizes are set correctly
    assert(n == (Index) model_->size_opt());
    assert(m == 0);

    // Get boundaries for the optimization parameters
    for(IdxT i = 0; i < model_->size_opt(); ++i)
    {
        x_l[i] = model_->param_lo()[i];
        x_u[i] = model_->param_up()[i];
    }

    return true;
}


template <class ScalarT, typename IdxT>
bool DynamicObjective<ScalarT, IdxT>::get_starting_point(Index n, bool init_x, Number* x,
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

    return true;
}


template <class ScalarT, typename IdxT>
bool DynamicObjective<ScalarT, IdxT>::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
    // Update optimization parameters
    for(IdxT i = 0; i < model_->size_opt(); ++i)
        model_->param()[i] = x[i];

    // Evaluate objective function
    integrator_->getSavedInitialCondition();
    integrator_->initializeSimulation(t_init_);
    integrator_->initializeQuadrature();
    integrator_->runSimulationQuadrature(t_final_, nout_);

    // Assuming objective function is given as the integral (quadrature) 0
    obj_value = (integrator_->getIntegral())[0];

    return true;
}


template <class ScalarT, typename IdxT>
bool DynamicObjective<ScalarT, IdxT>::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
    assert(model_->size_opt() == static_cast<IdxT>(n));
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
    integrator_->runForwardSimulation(t_final_, nout_);

    integrator_->initializeBackwardSimulation(t_final_);
    integrator_->runBackwardSimulation(t_init_);

    // For now assumes only one forward integrand and multiple optimization parameters.
    for(IdxT i = 0; i < model_->size_opt(); ++i)
        grad_f[i] = -((integrator_->getAdjointIntegral())[i]);

    integrator_->deleteAdjoint();

    return true;
}


template <class ScalarT, typename IdxT>
bool DynamicObjective<ScalarT, IdxT>::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
    return true;
}


template <class ScalarT, typename IdxT>
bool DynamicObjective<ScalarT, IdxT>::eval_jac_g(Index n, const Number* x, bool new_x,
                                                 Index m, Index nele_jac, Index* iRow, Index *jCol,
                                                 Number* values)
{
    return true;
}


template <class ScalarT, typename IdxT>
bool DynamicObjective<ScalarT, IdxT>::eval_h(Index n, const Number* x, bool new_x,
                                             Number obj_factor, Index m, const Number* lambda,
                                             bool new_lambda, Index nele_hess, Index* iRow,
                                             Index* jCol, Number* values)
{
    return true;
}


template <class ScalarT, typename IdxT>
void DynamicObjective<ScalarT, IdxT>::finalize_solution(SolverReturn status,
                                                        Index n, const Number* x, const Number* z_L, const Number* z_U,
                                                        Index m, const Number* g, const Number* lambda,
                                                        Number obj_value,
                                                        const IpoptData* ip_data,
                                                        IpoptCalculatedQuantities* ip_cq)
{

}




template class DynamicObjective<Ipopt::Number, Ipopt::Index>;
template class DynamicObjective<Ipopt::Number, std::size_t>;

} // namespace IpoptInterface
} // namespace AnalysisManager
