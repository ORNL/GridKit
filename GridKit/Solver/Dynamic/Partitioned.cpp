#include "Partitioned.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

#include <nvector/nvector_manyvector.h>
#include <sunlinsol/sunlinsol_spgmr.h>

#include <GridKit/Utilities/Logger/Logger.hpp>

#include <arkode/arkode.h>
#include <ida/ida.h>
#include <ida/ida_ls.h>

namespace AnalysisManager
{
  namespace Sundials
  {
    std::string PartitionedStats::report() const
    {
      static constexpr int label_width = 39;
      static constexpr int
                        stat_width = 12;
      std::stringstream out;
      out << std::setw(label_width) << "Outer steps" << " : "
          << std::setw(stat_width) << num_steps_ << '\n'
          << std::setw(label_width) << "Component residual evaluations" << " : "
          << std::setw(stat_width) << num_component_residual_evals_ << '\n'
          << std::setw(label_width) << "Component linear setups" << " : "
          << std::setw(stat_width) << num_component_linear_setups_ << '\n'
          << std::setw(label_width) << "Component error test failures" << " : "
          << std::setw(stat_width) << num_component_error_test_fails_ << '\n'
          << std::setw(label_width) << "Component nonlinear iterations" << " : "
          << std::setw(stat_width) << num_component_nonlinear_iters_ << '\n'
          << std::setw(label_width) << "Component nonlinear failures" << " : "
          << std::setw(stat_width) << num_component_nonlinear_failures_ << '\n'
          << std::setw(label_width) << "Coupling linear setups" << " : "
          << std::setw(stat_width) << num_coupling_linear_setups_ << '\n'
          << std::setw(label_width) << "Coupling nonlinear iterations" << " : "
          << std::setw(stat_width) << num_coupling_nonlinear_iters_ << '\n'
          << std::setw(label_width) << "Coupling nonlinear failures" << " : "
          << std::setw(stat_width) << num_coupling_nonlinear_failures_;
      return out.str();
    }

    template <class ScalarT, typename IdxT>
    Partitioned<ScalarT, IdxT>::Partitioned(EvaluatorT*       model,
                                            std::vector<Mask> component_masks,
                                            Mask              coupling_mask)
      : DynamicSolver<ScalarT, IdxT>(model),
        component_masks_(std::move(component_masks)),
        coupling_mask_(std::move(coupling_mask))
    {
      checkOutput(SUNContext_Create(SUN_COMM_NULL, &context_),
                  "SUNContext_Create");
    }

    template <class ScalarT, typename IdxT>
    Partitioned<ScalarT, IdxT>::~Partitioned()
    {
      deleteSimulation();
      SUNContext_Free(&context_);
    }

    template <class ScalarT, typename IdxT>
    int Partitioned<ScalarT, IdxT>::configureSimulation()
    {
      checkOutput(model_->initialize(), "Evaluator::initialize");
      checkOutput(model_->tagDifferentiable(), "Evaluator::tagDifferentiable");
      validateAndBuildIndices();
      allocateStateVectors();
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Partitioned<ScalarT, IdxT>::initializeSimulation(RealT t0)
    {
      // Honor state changes made between configureSimulation and initialization.
      for (size_t partition = 0; partition < component_masks_.size(); partition++)
      {
        gather(model_->y(), differential_indices_[partition], x_[partition]);
        gather(model_->y(), algebraic_indices_[partition], z_[partition]);
        gather(model_->yp(), differential_indices_[partition], xp_[partition]);
        gather(model_->yp(), algebraic_indices_[partition], zp_[partition]);
      }
      gather(model_->y(), coupling_indices_, w_);
      gather(model_->yp(), coupling_indices_, wp_);

      createSolver(t0);
      t_init_ = t0;

      if (model_->monitoring())
      {
        updateModelState(t0);
        model_->printMonitoredVariables();
      }
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Partitioned<ScalarT, IdxT>::runSimulation(
        RealT                                     tf,
        RealT                                     dt_monitor,
        std::optional<std::function<void(RealT)>> step_callback)
    {
      const int steps  = getMonitorStepCount(tf, dt_monitor);
      RealT     tret   = t_init_;
      int       retval = 0;
      for (int step = 1; step <= steps; step++)
      {
        const RealT tout = getMonitorTime(tf, dt_monitor, step, steps);
        retval           = ARKodeEvolve(solver_, tout, y_, &tret, ARK_NORMAL);
        checkOutput(retval, "ARKodeEvolve");

        if (step_callback.has_value() || model_->monitoring())
        {
          updateModelState(tret);
          if (model_->monitoring())
          {
            model_->printMonitoredVariables();
          }
          if (step_callback.has_value())
          {
            (*step_callback)(tret);
          }
        }
      }
      updateModelState(tret);
      return retval;
    }

    template <class ScalarT, typename IdxT>
    int Partitioned<ScalarT, IdxT>::deleteSimulation()
    {
      ARKodeFree(&solver_);
      SUNLinSolFree(algebraic_linear_solver_);
      algebraic_linear_solver_ = nullptr;
      for (auto& linear_solver : component_linear_solvers_)
      {
        SUNLinSolFree(linear_solver);
      }
      component_linear_solvers_.clear();

      N_VDestroy(y_);
      N_VDestroy(yp_);
      y_  = nullptr;
      yp_ = nullptr;
      for (auto& vector : x_)
      {
        N_VDestroy(vector);
      }
      for (auto& vector : z_)
      {
        N_VDestroy(vector);
      }
      for (auto& vector : xp_)
      {
        N_VDestroy(vector);
      }
      for (auto& vector : zp_)
      {
        N_VDestroy(vector);
      }
      N_VDestroy(w_);
      N_VDestroy(wp_);
      x_.clear();
      z_.clear();
      xp_.clear();
      zp_.clear();
      w_  = nullptr;
      wp_ = nullptr;
      return 0;
    }

    template <class ScalarT, typename IdxT>
    void Partitioned<ScalarT, IdxT>::setFixedStep(ScalarT time_step)
    {
      time_step_ = static_cast<RealT>(time_step);
    }

    template <class ScalarT, typename IdxT>
    void Partitioned<ScalarT, IdxT>::setTolerance(ScalarT rel_tol,
                                                  ScalarT abs_tol_override)
    {
      rel_tol_          = static_cast<RealT>(rel_tol);
      abs_tol_override_ = static_cast<RealT>(abs_tol_override);
    }

    template <class ScalarT, typename IdxT>
    void Partitioned<ScalarT, IdxT>::setMaxSteps(IdxT max_steps)
    {
      max_steps_ = max_steps;
    }

    template <class ScalarT, typename IdxT>
    void Partitioned<ScalarT, IdxT>::validateAndBuildIndices()
    {
      const size_t size = static_cast<size_t>(model_->size());
      if (coupling_mask_.size() != size)
      {
        throw PartitionedException("Partitioned coupling mask size does not match the model");
      }
      const auto& tags = model_->tag();
      if (tags.size() != size)
      {
        throw PartitionedException("Partitioned requires one differentiability tag per model variable");
      }

      // TODO: remove most of this when model supports partitioning
      differential_indices_.assign(component_masks_.size(), {});
      algebraic_indices_.assign(component_masks_.size(), {});
      component_indices_.assign(component_masks_.size(), {});
      coupling_indices_.clear();
      all_algebraic_indices_.clear();
      std::vector<int> owners(size, 0);

      for (size_t partition = 0; partition < component_masks_.size(); partition++)
      {
        const auto& mask = component_masks_[partition];
        if (mask.size() != size)
        {
          throw PartitionedException("Partitioned component mask size does not match the model");
        }
        for (size_t index = 0; index < size; index++)
        {
          if (!mask[index])
          {
            continue;
          }
          owners[index]++;
          const auto model_index = static_cast<IdxT>(index);
          if (tags[index])
          {
            differential_indices_[partition].push_back(model_index);
          }
          else
          {
            algebraic_indices_[partition].push_back(model_index);
          }
        }
        if (differential_indices_[partition].empty() && algebraic_indices_[partition].empty())
        {
          throw PartitionedException("Partitioned component masks may not be empty");
        }
        component_indices_[partition] = differential_indices_[partition];
        component_indices_[partition].insert(component_indices_[partition].end(),
                                             algebraic_indices_[partition].begin(),
                                             algebraic_indices_[partition].end());
        all_algebraic_indices_.insert(all_algebraic_indices_.end(),
                                      algebraic_indices_[partition].begin(),
                                      algebraic_indices_[partition].end());
      }

      for (size_t index = 0; index < size; index++)
      {
        if (coupling_mask_[index])
        {
          owners[index]++;
          if (tags[index])
          {
            throw PartitionedException("Partitioned coupling variables must be algebraic");
          }
          coupling_indices_.push_back(static_cast<IdxT>(index));
        }
        if (owners[index] != 1)
        {
          throw PartitionedException("Partitioned masks must cover every model index exactly once");
        }
      }
      if (coupling_indices_.empty())
      {
        throw PartitionedException("Partitioned requires at least one coupling variable");
      }
      all_algebraic_indices_.insert(all_algebraic_indices_.end(),
                                    coupling_indices_.begin(),
                                    coupling_indices_.end());
    }

    template <class ScalarT, typename IdxT>
    void Partitioned<ScalarT, IdxT>::allocateStateVectors()
    {
      const size_t partitions = component_masks_.size();
      x_.resize(partitions, nullptr);
      z_.resize(partitions, nullptr);
      xp_.resize(partitions, nullptr);
      zp_.resize(partitions, nullptr);
      for (size_t partition = 0; partition < partitions; partition++)
      {
        x_[partition] = N_VNew_Serial(
            static_cast<sunindextype>(differential_indices_[partition].size()), context_);
        z_[partition] = N_VNew_Serial(
            static_cast<sunindextype>(algebraic_indices_[partition].size()), context_);
        checkAllocation(x_[partition], "N_VNew_Serial");
        checkAllocation(z_[partition], "N_VNew_Serial");
        xp_[partition] = N_VClone(x_[partition]);
        zp_[partition] = N_VClone(z_[partition]);
        checkAllocation(xp_[partition], "N_VClone");
        checkAllocation(zp_[partition], "N_VClone");
      }
      w_ = N_VNew_Serial(static_cast<sunindextype>(coupling_indices_.size()), context_);
      checkAllocation(w_, "N_VNew_Serial");
      wp_ = N_VClone(w_);
      checkAllocation(wp_, "N_VClone");

      y_  = PDAEStepManyVector(x_.data(), z_.data(), w_, static_cast<int>(partitions));
      yp_ = PDAEStepManyVector(xp_.data(), zp_.data(), wp_, static_cast<int>(partitions));
      checkAllocation(y_, "PartitionedManyVector");
      checkAllocation(yp_, "PartitionedManyVector");
    }

    template <class ScalarT, typename IdxT>
    void Partitioned<ScalarT, IdxT>::createSolver(RealT t0)
    {
      solver_ = PDAEStepCreate(ComponentResidual, AlgebraicResidual, t0, y_, yp_, static_cast<int>(component_masks_.size()), context_);
      checkAllocation(solver_, "PartitionedCreate");
      checkOutput(ARKodeSetUserData(solver_, this), "ARKodeSetUserData");
      checkOutput(ARKodeSetFixedStep(solver_, time_step_), "ARKodeSetFixedStep");
      checkOutput(ARKodeSetMaxNumSteps(solver_, static_cast<long int>(max_steps_)),
                  "ARKodeSetMaxNumSteps");
      checkOutput(PDAEStepSetMaxNonlinIters(solver_, 16),
                  "PartitionedSetMaxNonlinIters");

      N_Vector algebraic_template = nullptr;
      checkOutput(PDAEStepGetAlgebraicVectorTemplate(solver_, &algebraic_template),
                  "PartitionedGetAlgebraicVectorTemplate");
      const int algebraic_maxl = static_cast<int>(N_VGetLength(algebraic_template));
      algebraic_linear_solver_ =
          SUNLinSol_SPGMR(algebraic_template, SUN_PREC_NONE, algebraic_maxl, context_);
      checkAllocation(algebraic_linear_solver_, "SUNLinSol_SPGMR");
      checkOutput(PDAEStepSetLinearSolver(solver_, algebraic_linear_solver_, nullptr),
                  "PartitionedSetLinearSolver");

      configurePartitionSolvers();
      configureTolerances();
      if (model_->hasJacobian())
      {
        checkOutput(PDAEStepSetPartitionJacTimes(solver_, ComponentJacTimes),
                    "PartitionedSetPartitionJacTimes");
        checkOutput(PDAEStepSetCouplingJacTimes(solver_, AlgebraicJacTimes),
                    "PartitionedSetCouplingJacTimes");
      }
    }

    template <class ScalarT, typename IdxT>
    void Partitioned<ScalarT, IdxT>::configurePartitionSolvers()
    {
      const size_t partitions = component_masks_.size();
      component_linear_solvers_.resize(partitions, nullptr);
      for (size_t partition = 0; partition < partitions; partition++)
      {
        N_Vector vector_template = nullptr;
        void*    ida_memory      = nullptr;
        checkOutput(PDAEStepGetPartitionVectorTemplate(
                        solver_, static_cast<int>(partition), &vector_template),
                    "PartitionedGetPartitionVectorTemplate");
        const int component_maxl = static_cast<int>(N_VGetLength(vector_template));
        component_linear_solvers_[partition] =
            SUNLinSol_SPGMR(vector_template, SUN_PREC_NONE, component_maxl, context_);
        checkAllocation(component_linear_solvers_[partition], "SUNLinSol_SPGMR");
        checkOutput(PDAEStepGetPartitionIntegrator(
                        solver_, static_cast<int>(partition), &ida_memory),
                    "PartitionedGetPartitionIntegrator");
        checkOutput(IDASetLinearSolver(ida_memory,
                                       component_linear_solvers_[partition],
                                       nullptr),
                    "IDASetLinearSolver");

        N_Vector id = N_VClone(vector_template);
        checkAllocation(id, "N_VClone");
        N_VConst(0.0, id);
        N_VConst(1.0, N_VGetSubvector_ManyVector(id, 0));
        checkOutput(IDASetId(ida_memory, id), "IDASetId");
        N_VDestroy(id);
      }
    }

    template <class ScalarT, typename IdxT>
    void Partitioned<ScalarT, IdxT>::configureTolerances()
    {
      if (abs_tol_override_ > 0)
      {
        checkOutput(ARKodeSStolerances(solver_, rel_tol_, abs_tol_override_),
                    "ARKodeSStolerances");
        for (size_t partition = 0; partition < component_masks_.size(); partition++)
        {
          void* ida_memory = nullptr;
          checkOutput(PDAEStepGetPartitionIntegrator(
                          solver_, static_cast<int>(partition), &ida_memory),
                      "PartitionedGetPartitionIntegrator");
          checkOutput(IDASStolerances(ida_memory, rel_tol_, abs_tol_override_),
                      "IDASStolerances");
        }
        return;
      }

      checkOutput(model_->setAbsoluteTolerance(rel_tol_),
                  "Evaluator::setAbsoluteTolerance");
      N_Vector absolute_tolerance = N_VClone(y_);
      checkAllocation(absolute_tolerance, "N_VClone");
      for (size_t partition = 0; partition < component_masks_.size(); partition++)
      {
        gather(model_->absoluteTolerance(), differential_indices_[partition], PDAEStepGetDifferentialSubvector(absolute_tolerance, static_cast<int>(partition)));
        gather(model_->absoluteTolerance(), algebraic_indices_[partition], PDAEStepGetAlgebraicSubvector(absolute_tolerance, static_cast<int>(partition)));
      }
      gather(model_->absoluteTolerance(), coupling_indices_, PDAEStepGetCouplingSubvector(absolute_tolerance));
      checkOutput(ARKodeSVtolerances(solver_, rel_tol_, absolute_tolerance),
                  "ARKodeSVtolerances");

      for (size_t partition = 0; partition < component_masks_.size(); partition++)
      {
        void*    ida_memory      = nullptr;
        N_Vector vector_template = nullptr;
        checkOutput(PDAEStepGetPartitionIntegrator(
                        solver_, static_cast<int>(partition), &ida_memory),
                    "PartitionedGetPartitionIntegrator");
        checkOutput(PDAEStepGetPartitionVectorTemplate(
                        solver_, static_cast<int>(partition), &vector_template),
                    "PartitionedGetPartitionVectorTemplate");
        N_Vector local_tolerance = N_VClone(vector_template);
        checkAllocation(local_tolerance, "N_VClone");
        gather(model_->absoluteTolerance(), differential_indices_[partition], N_VGetSubvector_ManyVector(local_tolerance, 0));
        gather(model_->absoluteTolerance(), algebraic_indices_[partition], N_VGetSubvector_ManyVector(local_tolerance, 1));
        checkOutput(IDASVtolerances(ida_memory, rel_tol_, local_tolerance),
                    "IDASVtolerances");
        N_VDestroy(local_tolerance);
      }
      N_VDestroy(absolute_tolerance);
    }

    template <class ScalarT, typename IdxT>
    int Partitioned<ScalarT, IdxT>::ComponentResidual(
        int partition, RealT t, N_Vector y, N_Vector w, N_Vector yp, N_Vector residual, void* user_data)
    {
      try
      {
        return static_cast<Partitioned*>(user_data)->componentResidual(
            partition, t, y, w, yp, residual);
      }
      catch (...)
      {
        return -1;
      }
    }

    template <class ScalarT, typename IdxT>
    int Partitioned<ScalarT, IdxT>::ComponentJacTimes(
        int partition, RealT t, N_Vector y, N_Vector w, N_Vector yp, N_Vector, N_Vector v, N_Vector Jv, RealT cj, void* user_data, N_Vector, N_Vector)
    {
      try
      {
        return static_cast<Partitioned*>(user_data)->componentJacTimes(
            partition, t, cj, y, w, yp, v, Jv);
      }
      catch (...)
      {
        return -1;
      }
    }

    template <class ScalarT, typename IdxT>
    int Partitioned<ScalarT, IdxT>::AlgebraicResidual(
        RealT t, N_Vector y, N_Vector w, N_Vector residual, void* user_data)
    {
      try
      {
        return static_cast<Partitioned*>(user_data)->algebraicResidual(t, y, w, residual);
      }
      catch (...)
      {
        return -1;
      }
    }

    template <class ScalarT, typename IdxT>
    int Partitioned<ScalarT, IdxT>::AlgebraicJacTimes(
        RealT t, N_Vector y, N_Vector w, N_Vector v, N_Vector Jv, void* user_data, N_Vector)
    {
      try
      {
        return static_cast<Partitioned*>(user_data)->algebraicJacTimes(t, y, w, v, Jv);
      }
      catch (...)
      {
        return -1;
      }
    }

    template <class ScalarT, typename IdxT>
    int Partitioned<ScalarT, IdxT>::componentResidual(
        int partition, RealT t, N_Vector y, N_Vector w, N_Vector yp, N_Vector residual)
    {
      const auto part = static_cast<size_t>(partition);
      // TODO: remove scattering
      scatterComponentState(partition, y, w, yp);
      model_->updateTime(t, 0.0);
      const int retval = model_->evaluateResidual();
      if (retval != 0)
      {
        return retval;
      }
      gather(model_->getResidual(), differential_indices_[part], N_VGetSubvector_ManyVector(residual, 0));
      gather(model_->getResidual(), algebraic_indices_[part], N_VGetSubvector_ManyVector(residual, 1));
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Partitioned<ScalarT, IdxT>::componentJacTimes(
        int partition, RealT t, RealT cj, N_Vector y, N_Vector w, N_Vector yp, N_Vector v, N_Vector Jv)
    {
      const auto part = static_cast<size_t>(partition);
      // TODO: remove scattering
      scatterComponentState(partition, y, w, yp);
      model_->updateTime(t, cj);
      const int retval = model_->evaluateJacobian();
      if (retval != 0)
      {
        return retval;
      }
      multiplyMaskedJacobian(component_indices_[part], v, Jv);
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Partitioned<ScalarT, IdxT>::algebraicResidual(
        RealT t, N_Vector y, N_Vector w, N_Vector residual)
    {
      // TODO: remove scattering
      scatterAlgebraicState(y, w);
      model_->updateTime(t, 0.0);
      const int retval = model_->evaluateResidual();
      if (retval != 0)
      {
        return retval;
      }
      for (size_t partition = 0; partition < algebraic_indices_.size(); partition++)
      {
        auto subvector = N_VGetSubvector_ManyVector(
            residual, static_cast<sunindextype>(partition));
        gather(model_->getResidual(), algebraic_indices_[partition], subvector);
      }
      gather(model_->getResidual(), coupling_indices_, N_VGetSubvector_ManyVector(residual, static_cast<sunindextype>(algebraic_indices_.size())));
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Partitioned<ScalarT, IdxT>::algebraicJacTimes(
        RealT t, N_Vector y, N_Vector w, N_Vector v, N_Vector Jv)
    {
      // TODO: remove scattering
      scatterAlgebraicState(y, w);
      model_->updateTime(t, 0.0);
      const int retval = model_->evaluateJacobian();
      if (retval != 0)
      {
        return retval;
      }
      multiplyMaskedJacobian(all_algebraic_indices_, v, Jv);
      return 0;
    }

    template <class ScalarT, typename IdxT>
    void Partitioned<ScalarT, IdxT>::scatterComponentState(
        int partition, N_Vector y, N_Vector w, N_Vector yp)
    {
      const auto part = static_cast<size_t>(partition);
      scatter(N_VGetSubvector_ManyVector(y, 0), differential_indices_[part], model_->y());
      scatter(N_VGetSubvector_ManyVector(y, 1), algebraic_indices_[part], model_->y());
      scatter(w, coupling_indices_, model_->y());
      scatter(N_VGetSubvector_ManyVector(yp, 0), differential_indices_[part], model_->yp());
      scatter(N_VGetSubvector_ManyVector(yp, 1), algebraic_indices_[part], model_->yp());
    }

    template <class ScalarT, typename IdxT>
    void Partitioned<ScalarT, IdxT>::scatterFullState(N_Vector y)
    {
      // TODO: remove scattering
      for (size_t partition = 0; partition < component_masks_.size(); partition++)
      {
        scatter(PDAEStepGetDifferentialSubvector(y, static_cast<int>(partition)),
                differential_indices_[partition],
                model_->y());
        scatter(PDAEStepGetAlgebraicSubvector(y, static_cast<int>(partition)),
                algebraic_indices_[partition],
                model_->y());
      }
      scatter(PDAEStepGetCouplingSubvector(y), coupling_indices_, model_->y());
    }

    template <class ScalarT, typename IdxT>
    void Partitioned<ScalarT, IdxT>::scatterAlgebraicState(N_Vector y, N_Vector w)
    {
      // TODO: remove scattering
      for (size_t partition = 0; partition < component_masks_.size(); partition++)
      {
        scatter(PDAEStepGetDifferentialSubvector(y, static_cast<int>(partition)),
                differential_indices_[partition],
                model_->y());
        scatter(PDAEStepGetAlgebraicSubvector(y, static_cast<int>(partition)),
                algebraic_indices_[partition],
                model_->y());
      }
      scatter(w, coupling_indices_, model_->y());
    }

    template <class ScalarT, typename IdxT>
    void Partitioned<ScalarT, IdxT>::multiplyMaskedJacobian(
        const std::vector<IdxT>& indices, N_Vector v, N_Vector Jv)
    {
      // TODO: remove masking
      std::vector<RealT> global_v(static_cast<size_t>(model_->size()), 0.0);
      size_t             position   = 0;
      const auto         subvectors = N_VGetNumSubvectors_ManyVector(v);
      for (sunindextype sub = 0; sub < subvectors; sub++)
      {
        N_Vector     serial = N_VGetSubvector_ManyVector(v, sub);
        const auto   length = static_cast<size_t>(N_VGetLength(serial));
        const RealT* data   = N_VGetArrayPointer(serial);
        for (size_t local = 0; local < length; local++)
        {
          global_v[static_cast<size_t>(indices[position++])] = data[local];
        }
      }

      auto* jacobian = model_->getCsrJacobian();
      if (jacobian == nullptr)
      {
        throw PartitionedException("Evaluator reported a Jacobian but returned null");
      }
      const IdxT*        rows = jacobian->getRowData();
      const IdxT*        cols = jacobian->getColData();
      const RealT*       vals = jacobian->getValues();
      std::vector<RealT> result(indices.size(), 0.0);
      for (size_t local_row = 0; local_row < indices.size(); local_row++)
      {
        const IdxT row = indices[local_row];
        for (IdxT entry = rows[row]; entry < rows[row + 1]; entry++)
        {
          result[local_row] += vals[entry] * global_v[static_cast<size_t>(cols[entry])];
        }
      }

      position                     = 0;
      const auto output_subvectors = N_VGetNumSubvectors_ManyVector(Jv);
      for (sunindextype sub = 0; sub < output_subvectors; sub++)
      {
        N_Vector   serial = N_VGetSubvector_ManyVector(Jv, sub);
        const auto length = static_cast<size_t>(N_VGetLength(serial));
        RealT*     data   = N_VGetArrayPointer(serial);
        for (size_t local = 0; local < length; local++)
        {
          data[local] = result[position++];
        }
      }
    }

    template <class ScalarT, typename IdxT>
    void Partitioned<ScalarT, IdxT>::updateModelState(RealT t)
    {
      // TODO: significant updates once the model supports partitioning
      scatterFullState(y_);
      for (size_t partition = 0; partition < component_masks_.size(); partition++)
      {
        void*    ida_memory = nullptr;
        N_Vector current_yp = nullptr;
        checkOutput(PDAEStepGetPartitionIntegrator(
                        solver_, static_cast<int>(partition), &ida_memory),
                    "PartitionedGetPartitionIntegrator");
        checkOutput(IDAGetCurrentYp(ida_memory, &current_yp), "IDAGetCurrentYp");
        // IDA does not publish its current derivative vector until its first
        // evolve call.  Monitoring at the initial time must therefore use the
        // derivative ManyVector supplied to PDAEStepCreate.
        if (current_yp == nullptr)
        {
          scatter(PDAEStepGetDifferentialSubvector(
                      yp_, static_cast<int>(partition)),
                  differential_indices_[partition],
                  model_->yp());
          scatter(PDAEStepGetAlgebraicSubvector(
                      yp_, static_cast<int>(partition)),
                  algebraic_indices_[partition],
                  model_->yp());
        }
        else
        {
          scatter(N_VGetSubvector_ManyVector(current_yp, 0),
                  differential_indices_[partition],
                  model_->yp());
          scatter(N_VGetSubvector_ManyVector(current_yp, 1),
                  algebraic_indices_[partition],
                  model_->yp());
        }
      }
      auto* derivative = model_->yp().getData();
      for (const IdxT index : coupling_indices_)
      {
        derivative[index] = 0.0;
      }
      model_->yp().setDataUpdated();
      model_->updateTime(t, 0.0);
    }

    template <class ScalarT, typename IdxT>
    void Partitioned<ScalarT, IdxT>::gather(
        const VectorT& source, const std::vector<IdxT>& indices, N_Vector destination) const
    {
      const ScalarT* source_data      = source.getData();
      RealT*         destination_data = N_VGetArrayPointer(destination);
      for (size_t local = 0; local < indices.size(); local++)
      {
        destination_data[local] = static_cast<RealT>(source_data[indices[local]]);
      }
    }

    template <class ScalarT, typename IdxT>
    void Partitioned<ScalarT, IdxT>::scatter(
        N_Vector source, const std::vector<IdxT>& indices, VectorT& destination) const
    {
      const RealT* source_data      = N_VGetArrayPointer(source);
      ScalarT*     destination_data = destination.getData();
      for (size_t local = 0; local < indices.size(); local++)
      {
        destination_data[indices[local]] = static_cast<ScalarT>(source_data[local]);
      }
      destination.setDataUpdated();
    }

    template <class ScalarT, typename IdxT>
    PartitionedStats Partitioned<ScalarT, IdxT>::getStats() const
    {
      PartitionedStats stats;
      checkOutput(ARKodeGetNumSteps(solver_, &stats.num_steps_),
                  "ARKodeGetNumSteps");
      checkOutput(PDAEStepGetNumLinSolvSetups(
                      solver_, &stats.num_coupling_linear_setups_),
                  "PartitionedGetNumLinSolvSetups");
      checkOutput(PDAEStepGetNonlinSolvStats(
                      solver_, &stats.num_coupling_nonlinear_iters_, &stats.num_coupling_nonlinear_failures_),
                  "PartitionedGetNonlinSolvStats");
      for (size_t partition = 0; partition < component_masks_.size(); partition++)
      {
        void*    ida_memory           = nullptr;
        long int nonlinear_iterations = 0;
        long int nonlinear_failures   = 0;
        long int value                = 0;
        checkOutput(PDAEStepGetPartitionIntegrator(
                        solver_, static_cast<int>(partition), &ida_memory),
                    "PartitionedGetPartitionIntegrator");
        checkOutput(IDAGetNumResEvals(ida_memory, &value), "IDAGetNumResEvals");
        stats.num_component_residual_evals_ += value;
        checkOutput(IDAGetNumLinSolvSetups(ida_memory, &value),
                    "IDAGetNumLinSolvSetups");
        stats.num_component_linear_setups_ += value;
        checkOutput(IDAGetNumErrTestFails(ida_memory, &value),
                    "IDAGetNumErrTestFails");
        stats.num_component_error_test_fails_ += value;
        checkOutput(IDAGetNonlinSolvStats(ida_memory, &nonlinear_iterations, &nonlinear_failures),
                    "IDAGetNonlinSolvStats");
        stats.num_component_nonlinear_iters_    += nonlinear_iterations;
        stats.num_component_nonlinear_failures_ += nonlinear_failures;
      }
      return stats;
    }

    template <class ScalarT, typename IdxT>
    int Partitioned<ScalarT, IdxT>::getMonitorStepCount(RealT tf,
                                                        RealT dt_monitor) const
    {
      if (dt_monitor <= 0.0)
      {
        return 1;
      }
      const RealT estimate = (tf - t_init_) / dt_monitor;
      const RealT epsilon  = std::numeric_limits<RealT>::epsilon()
                            * std::max({std::abs(t_init_), std::abs(tf), RealT(1.0)})
                            / dt_monitor;
      return static_cast<int>(std::ceil(estimate - epsilon));
    }

    template <class ScalarT, typename IdxT>
    typename Partitioned<ScalarT, IdxT>::RealT
    Partitioned<ScalarT, IdxT>::getMonitorTime(RealT tf, RealT dt_monitor, int step, int steps) const
    {
      return step == steps ? tf : std::fma(static_cast<RealT>(step), dt_monitor, t_init_);
    }

    template <class ScalarT, typename IdxT>
    void Partitioned<ScalarT, IdxT>::checkOutput(int         retval,
                                                 const char* function_name)
    {
      if (retval < 0)
      {
        GridKit::Utilities::Logger::error()
            << function_name << " failed with flag " << retval << ".\n";
        throw PartitionedException(std::string(function_name) + " failed");
      }
    }

    template <class ScalarT, typename IdxT>
    void Partitioned<ScalarT, IdxT>::checkAllocation(
        const void* pointer, const char* function_name)
    {
      if (pointer == nullptr)
      {
        GridKit::Utilities::Logger::error()
            << function_name << " returned a null pointer.\n";
        throw PartitionedException(std::string(function_name) + " failed");
      }
    }

    template class Partitioned<sunrealtype, long int>;
    template class Partitioned<sunrealtype, int>;
    template class Partitioned<sunrealtype, size_t>;
  } // namespace Sundials
} // namespace AnalysisManager
