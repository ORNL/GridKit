
#pragma once

#include <cassert>
#include <iostream>
#include <vector>

#include <GridKit/Model/PowerFlow/ModelEvaluatorImpl.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{

  /**
   * @brief Prototype for a system model class
   *
   * This class maps component data to system data and implements
   * Model::Evaluator for the system model. This is still work in
   * progress and code is not optimized.
   *
   * @todo Address thread safety for the system model methods.
   *
   */
  template <class ScalarT, typename IdxT>
  class SystemModel : public ModelEvaluatorImpl<ScalarT, IdxT>
  {
    using bus_type       = Model::Evaluator<ScalarT, IdxT>;
    using component_type = Model::Evaluator<ScalarT, IdxT>;
    using RealT          = typename ModelEvaluatorImpl<ScalarT, IdxT>::RealT;
    using VectorT        = typename ModelEvaluatorImpl<ScalarT, IdxT>::VectorT;

    using ModelEvaluatorImpl<ScalarT, IdxT>::size_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::size_quad_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::size_opt_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::nnz_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::time_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::alpha_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::y_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::yp_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::yB_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::ypB_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::tag_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::abs_tol_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::f_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::fB_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::g_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::gB_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::param_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::param_up_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::param_lo_;

  public:
    /**
     * @brief Constructor for the system model
     */
    SystemModel()
      : ModelEvaluatorImpl<ScalarT, IdxT>(0, 0, 0)
    {
    }

    /**
     * @brief Destructor for the system model
     */
    virtual ~SystemModel()
    {
    }

    /**
     * @brief Allocate buses, components, and system objects.
     *
     * This method first allocates bus objects, then component objects,
     * and computes system size (number of unknowns). Once the size is
     * computed, system global objects are allocated.
     *
     * @post size_quad_ == 0 or 1
     * @post size_ >= 1
     * @post size_opt_ >= 0
     *
     */
    int allocate()
    {
      size_      = 0;
      size_quad_ = 0;
      size_opt_  = 0;

      // Allocate all buses
      for (const auto& bus : buses_)
      {
        bus->allocate();
        size_      += bus->size();
        size_quad_ += bus->sizeQuadrature();
        size_opt_  += bus->sizeParams();
      }

      // Allocate all components
      for (const auto& component : components_)
      {
        component->allocate();
        size_      += component->size();
        size_quad_ += component->sizeQuadrature();
        size_opt_  += component->sizeParams();
      }

      // Allocate global vectors
      this->allocateVectors(size_);
      tag_.resize(size_);

      yB_.resize(size_);
      ypB_.resize(size_);
      fB_.resize(size_);
      g_.resize(size_quad_);
      gB_.resize(size_quad_ * size_opt_);
      param_.resize(size_opt_);
      param_lo_.resize(size_opt_);
      param_up_.resize(size_opt_);

      assert(size_quad_ == 1 or size_quad_ == 0);

      return 0;
    }

    /**
     * @brief Assume that jacobian is not avalible
     *
     * @return true
     * @return false
     */
    bool hasJacobian()
    {
      return false;
    }

    /**
     * @brief Initialize buses first, then all the other components.
     *
     * @pre All buses and components must be allocated at this point.
     * @pre Bus variables are written before component variables in the
     * system variable vector.
     *
     * Buses must be initialized before other components, because other
     * components may write to buses during the initialization.
     *
     * Also, generators may write to control devices (e.g. governors,
     * exciters, etc.) during the initialization.
     *
     * @todo Implement writting to system vectors in a thread-safe way.
     */
    int initialize()
    {
      // Set initial values for global solution vectors
      IdxT     varOffset = 0;
      IdxT     optOffset = 0;
      auto*    y         = y_.getData();
      auto*    yp        = yp_.getData();
      ScalarT* param     = nullptr;
      ScalarT* param_lo  = nullptr;
      ScalarT* param_up  = nullptr;

      if (this->sizeParams() > 0)
      {
        param    = param_.getData();
        param_lo = param_lo_.getData();
        param_up = param_up_.getData();
      }

      for (const auto& bus : buses_)
      {
        bus->initialize();
      }

      for (const auto& bus : buses_)
      {
        if (bus->size() > 0)
        {
          const auto* bus_y  = bus->y().getData();
          const auto* bus_yp = bus->yp().getData();

          for (IdxT j = 0; j < bus->size(); ++j)
          {
            y[varOffset + j]  = bus_y[j];
            yp[varOffset + j] = bus_yp[j];
          }
        }
        varOffset += bus->size();

        if (bus->sizeParams() > 0)
        {
          const auto* bus_param    = bus->param().getData();
          const auto* bus_param_lo = bus->param_lo().getData();
          const auto* bus_param_up = bus->param_up().getData();

          for (IdxT j = 0; j < bus->sizeParams(); ++j)
          {
            param[optOffset + j]    = bus_param[j];
            param_lo[optOffset + j] = bus_param_lo[j];
            param_up[optOffset + j] = bus_param_up[j];
          }
        }
        optOffset += bus->sizeParams();
      }

      // Initialize components
      for (const auto& component : components_)
      {
        component->initialize();
      }

      for (const auto& component : components_)
      {
        if (component->size() > 0)
        {
          const auto* component_y  = component->y().getData();
          const auto* component_yp = component->yp().getData();

          for (IdxT j = 0; j < component->size(); ++j)
          {
            y[varOffset + j]  = component_y[j];
            yp[varOffset + j] = component_yp[j];
          }
        }
        varOffset += component->size();

        if (component->sizeParams() > 0)
        {
          const auto* component_param    = component->param().getData();
          const auto* component_param_lo = component->param_lo().getData();
          const auto* component_param_up = component->param_up().getData();

          for (IdxT j = 0; j < component->sizeParams(); ++j)
          {
            param[optOffset + j]    = component_param[j];
            param_lo[optOffset + j] = component_param_lo[j];
            param_up[optOffset + j] = component_param_up[j];
          }
        }
        optOffset += component->sizeParams();
      }

      y_.setDataUpdated();
      yp_.setDataUpdated();
      if (this->sizeParams() > 0)
      {
        param_.setDataUpdated();
        param_lo_.setDataUpdated();
        param_up_.setDataUpdated();
      }

      return 0;
    }

    /**
     * @todo Tagging differential variables
     *
     * Identify what variables in the system of differential-algebraic
     * equations are differential variables, i.e. their derivatives
     * appear in the equations.
     */
    int tagDifferentiable()
    {
      // Set initial values for global solution vectors
      IdxT offset = 0;

      for (const auto& bus : buses_)
      {
        bus->tagDifferentiable();
        for (IdxT j = 0; j < bus->size(); ++j)
        {
          tag_[offset + j] = bus->tag()[j];
        }
        offset += bus->size();
      }

      for (const auto& component : components_)
      {
        component->tagDifferentiable();
        for (IdxT j = 0; j < component->size(); ++j)
        {
          tag_[offset + j] = component->tag()[j];
        }
        offset += component->size();
      }

      abs_tol_.setDataUpdated();

      return 0;
    }

    /**
     * @brief Compute the absolute tolerance for each variable in the model
     *
     * @param rel_tol The relative tolerance which can be used to pick the
     *        absolute tolerance.
     * @tparam ScalarT Scalar data type
     * @tparam IdxT Index data type
     * @return int 0 if successful, non-zero otherwise.
     *
     * This represents a "noise" level close to zero for which pure relative
     * error cannot be used.
     */
    int setAbsoluteTolerance(RealT rel_tol)
    {
      // Set initial values for global solution vectors
      IdxT  offset  = 0;
      auto* abs_tol = abs_tol_.getData();

      for (const auto& bus : buses_)
      {
        bus->setAbsoluteTolerance(rel_tol);
        if (bus->size() > 0)
        {
          const auto* bus_abs_tol = bus->absoluteTolerance().getData();
          for (IdxT j = 0; j < bus->size(); ++j)
          {
            abs_tol[offset + j] = bus_abs_tol[j];
          }
        }
        offset += bus->size();
      }

      for (const auto& component : components_)
      {
        component->setAbsoluteTolerance(rel_tol);
        if (component->size() > 0)
        {
          const auto* component_abs_tol = component->absoluteTolerance().getData();
          for (IdxT j = 0; j < component->size(); ++j)
          {
            abs_tol[offset + j] = component_abs_tol[j];
          }
        }
        offset += component->size();
      }

      abs_tol_.setDataUpdated();

      return 0;
    }

    /**
     * @brief Compute system residual vector
     *
     * First, update bus and component variables from the system solution
     * vector. Next, evaluate residuals in buses and components, and
     * then copy values to the global residual vector.
     *
     * @warning Residuals must be computed for buses, before component
     * residuals are computed. Buses own residuals for active and
     * power P and Q, but the contributions to these residuals come
     * from components. Buses assign their residual values, while components
     * add to those values by in-place adition. This is why bus residuals
     * need to be computed first.
     *
     * @todo Here, components write to local values, which are then copied
     * to global system vectors. Make components write to the system
     * vectors directly.
     */
    int evaluateResidual()
    {
      // Update variables
      IdxT           varOffset = 0;
      IdxT           optOffset = 0;
      const auto*    y         = y_.getData();
      const auto*    yp        = yp_.getData();
      auto*          f         = f_.getData();
      const ScalarT* param     = nullptr;

      if (this->sizeParams() > 0)
      {
        param = param_.getData();
      }

      for (const auto& bus : buses_)
      {
        if (bus->size() > 0)
        {
          auto* bus_y  = bus->y().getData();
          auto* bus_yp = bus->yp().getData();

          for (IdxT j = 0; j < bus->size(); ++j)
          {
            bus_y[j]  = y[varOffset + j];
            bus_yp[j] = yp[varOffset + j];
          }
          bus->y().setDataUpdated();
          bus->yp().setDataUpdated();
        }
        varOffset += bus->size();

        if (bus->sizeParams() > 0)
        {
          auto* bus_param = bus->param().getData();
          for (IdxT j = 0; j < bus->sizeParams(); ++j)
          {
            bus_param[j] = param[optOffset + j];
          }
          bus->param().setDataUpdated();
        }
        optOffset += bus->sizeParams();

        bus->evaluateResidual();
      }

      for (const auto& component : components_)
      {
        if (component->size() > 0)
        {
          auto* component_y  = component->y().getData();
          auto* component_yp = component->yp().getData();

          for (IdxT j = 0; j < component->size(); ++j)
          {
            component_y[j]  = y[varOffset + j];
            component_yp[j] = yp[varOffset + j];
          }
          component->y().setDataUpdated();
          component->yp().setDataUpdated();
        }
        varOffset += component->size();

        if (component->sizeParams() > 0)
        {
          auto* component_param = component->param().getData();
          for (IdxT j = 0; j < component->sizeParams(); ++j)
          {
            component_param[j] = param[optOffset + j];
          }
          component->param().setDataUpdated();
        }
        optOffset += component->sizeParams();

        component->evaluateResidual();
      }

      // Update residual vector
      IdxT resOffset = 0;
      for (const auto& bus : buses_)
      {
        if (bus->size() > 0)
        {
          const auto* bus_f = bus->getResidual().getData();
          for (IdxT j = 0; j < bus->size(); ++j)
          {
            f[resOffset + j] = bus_f[j];
          }
        }
        resOffset += bus->size();
      }

      for (const auto& component : components_)
      {
        if (component->size() > 0)
        {
          const auto* component_f = component->getResidual().getData();
          for (IdxT j = 0; j < component->size(); ++j)
          {
            f[resOffset + j] = component_f[j];
          }
        }
        resOffset += component->size();
      }

      f_.setDataUpdated();

      return 0;
    }

    /**
     * @brief Evaluate system Jacobian.
     *
     * @todo Need to implement Jacobian. For now, using finite difference
     * approximation provided by IDA. This works for dense Jacobian matrix
     * only.
     *
     */
    int evaluateJacobian()
    {
      return 0;
    }

    /**
     * @brief Evaluate integrands for the system quadratures.
     */
    int evaluateIntegrand()
    {
      // Update variables
      IdxT           varOffset = 0;
      IdxT           optOffset = 0;
      const auto*    y         = y_.getData();
      const auto*    yp        = yp_.getData();
      const ScalarT* param     = nullptr;
      ScalarT*       g         = nullptr;

      if (this->sizeParams() > 0)
      {
        param = param_.getData();
      }
      if (this->sizeQuadrature() > 0)
      {
        g = g_.getData();
      }

      for (const auto& bus : buses_)
      {
        if (bus->size() > 0)
        {
          auto* bus_y  = bus->y().getData();
          auto* bus_yp = bus->yp().getData();

          for (IdxT j = 0; j < bus->size(); ++j)
          {
            bus_y[j]  = y[varOffset + j];
            bus_yp[j] = yp[varOffset + j];
          }
          bus->y().setDataUpdated();
          bus->yp().setDataUpdated();
        }
        varOffset += bus->size();

        if (bus->sizeParams() > 0)
        {
          auto* bus_param = bus->param().getData();
          for (IdxT j = 0; j < bus->sizeParams(); ++j)
          {
            bus_param[j] = param[optOffset + j];
          }
          bus->param().setDataUpdated();
        }
        optOffset += bus->sizeParams();

        bus->evaluateIntegrand();
      }

      for (const auto& component : components_)
      {
        if (component->size() > 0)
        {
          auto* component_y  = component->y().getData();
          auto* component_yp = component->yp().getData();

          for (IdxT j = 0; j < component->size(); ++j)
          {
            component_y[j]  = y[varOffset + j];
            component_yp[j] = yp[varOffset + j];
          }
          component->y().setDataUpdated();
          component->yp().setDataUpdated();
        }
        varOffset += component->size();

        if (component->sizeParams() > 0)
        {
          auto* component_param = component->param().getData();
          for (IdxT j = 0; j < component->sizeParams(); ++j)
          {
            component_param[j] = param[optOffset + j];
          }
          component->param().setDataUpdated();
        }
        optOffset += component->sizeParams();

        component->evaluateIntegrand();
      }

      // Update integrand vector
      IdxT intOffset = 0;
      for (const auto& bus : buses_)
      {
        if (bus->sizeQuadrature() > 0)
        {
          const auto* bus_g = bus->getIntegrand().getData();
          for (IdxT j = 0; j < bus->sizeQuadrature(); ++j)
          {
            g[intOffset + j] = bus_g[j];
          }
        }
        intOffset += bus->sizeQuadrature();
      }

      for (const auto& component : components_)
      {
        if (component->sizeQuadrature() > 0)
        {
          const auto* component_g = component->getIntegrand().getData();
          for (IdxT j = 0; j < component->sizeQuadrature(); ++j)
          {
            g[intOffset + j] = component_g[j];
          }
        }
        intOffset += component->sizeQuadrature();
      }

      if (this->sizeQuadrature() > 0)
      {
        g_.setDataUpdated();
      }

      return 0;
    }

    /**
     * @brief Initialize system adjoint.
     *
     * Updates variables and optimization parameters, then initializes
     * adjoints locally and copies them to the system adjoint vector.
     */
    int initializeAdjoint()
    {
      IdxT           offset    = 0;
      IdxT           optOffset = 0;
      const auto*    y         = y_.getData();
      const auto*    yp        = yp_.getData();
      auto*          yB        = yB_.getData();
      auto*          ypB       = ypB_.getData();
      const ScalarT* param     = nullptr;

      if (this->sizeParams() > 0)
      {
        param = param_.getData();
      }

      // Update bus variables and optimization parameters
      for (const auto& bus : buses_)
      {
        if (bus->size() > 0)
        {
          auto* bus_y  = bus->y().getData();
          auto* bus_yp = bus->yp().getData();

          for (IdxT j = 0; j < bus->size(); ++j)
          {
            bus_y[j]  = y[offset + j];
            bus_yp[j] = yp[offset + j];
          }
          bus->y().setDataUpdated();
          bus->yp().setDataUpdated();
        }
        offset += bus->size();

        if (bus->sizeParams() > 0)
        {
          auto* bus_param = bus->param().getData();
          for (IdxT j = 0; j < bus->sizeParams(); ++j)
          {
            bus_param[j] = param[optOffset + j];
          }
          bus->param().setDataUpdated();
        }
        optOffset += bus->sizeParams();
      }

      // Update component variables and optimization parameters
      for (const auto& component : components_)
      {
        if (component->size() > 0)
        {
          auto* component_y  = component->y().getData();
          auto* component_yp = component->yp().getData();

          for (IdxT j = 0; j < component->size(); ++j)
          {
            component_y[j]  = y[offset + j];
            component_yp[j] = yp[offset + j];
          }
          component->y().setDataUpdated();
          component->yp().setDataUpdated();
        }
        offset += component->size();

        if (component->sizeParams() > 0)
        {
          auto* component_param = component->param().getData();
          for (IdxT j = 0; j < component->sizeParams(); ++j)
          {
            component_param[j] = param[optOffset + j];
          }
          component->param().setDataUpdated();
        }
        optOffset += component->sizeParams();
      }

      // Reset counter
      offset = 0;

      // Initialize bus adjoints
      for (const auto& bus : buses_)
      {
        bus->initializeAdjoint();
        if (bus->size() > 0)
        {
          const auto* bus_yB  = bus->yB().getData();
          const auto* bus_ypB = bus->ypB().getData();

          for (IdxT j = 0; j < bus->size(); ++j)
          {
            yB[offset + j]  = bus_yB[j];
            ypB[offset + j] = bus_ypB[j];
          }
        }
        offset += bus->size();
      }

      // Initialize component adjoints
      for (const auto& component : components_)
      {
        component->initializeAdjoint();
        if (component->size() > 0)
        {
          const auto* component_yB  = component->yB().getData();
          const auto* component_ypB = component->ypB().getData();

          for (IdxT j = 0; j < component->size(); ++j)
          {
            yB[offset + j]  = component_yB[j];
            ypB[offset + j] = component_ypB[j];
          }
        }
        offset += component->size();
      }

      yB_.setDataUpdated();
      ypB_.setDataUpdated();

      return 0;
    }

    /**
     * @brief Compute adjoint residual for the system model.
     *
     * @warning Components write to bus residuals. Do not copy bus residuals
     * to system vectors before components computed their residuals.
     *
     */
    int evaluateAdjointResidual()
    {
      IdxT           varOffset = 0;
      IdxT           optOffset = 0;
      const auto*    y         = y_.getData();
      const auto*    yp        = yp_.getData();
      const auto*    yB        = yB_.getData();
      const auto*    ypB       = ypB_.getData();
      auto*          fB        = fB_.getData();
      const ScalarT* param     = nullptr;

      if (this->sizeParams() > 0)
      {
        param = param_.getData();
      }

      // Update variables in component models
      for (const auto& bus : buses_)
      {
        if (bus->size() > 0)
        {
          auto* bus_y   = bus->y().getData();
          auto* bus_yp  = bus->yp().getData();
          auto* bus_yB  = bus->yB().getData();
          auto* bus_ypB = bus->ypB().getData();

          for (IdxT j = 0; j < bus->size(); ++j)
          {
            bus_y[j]   = y[varOffset + j];
            bus_yp[j]  = yp[varOffset + j];
            bus_yB[j]  = yB[varOffset + j];
            bus_ypB[j] = ypB[varOffset + j];
          }
          bus->y().setDataUpdated();
          bus->yp().setDataUpdated();
          bus->yB().setDataUpdated();
          bus->ypB().setDataUpdated();
        }
        varOffset += bus->size();

        if (bus->sizeParams() > 0)
        {
          auto* bus_param = bus->param().getData();
          for (IdxT j = 0; j < bus->sizeParams(); ++j)
          {
            bus_param[j] = param[optOffset + j];
          }
          bus->param().setDataUpdated();
        }
        optOffset += bus->sizeParams();
      }

      for (const auto& component : components_)
      {
        if (component->size() > 0)
        {
          auto* component_y   = component->y().getData();
          auto* component_yp  = component->yp().getData();
          auto* component_yB  = component->yB().getData();
          auto* component_ypB = component->ypB().getData();

          for (IdxT j = 0; j < component->size(); ++j)
          {
            component_y[j]   = y[varOffset + j];
            component_yp[j]  = yp[varOffset + j];
            component_yB[j]  = yB[varOffset + j];
            component_ypB[j] = ypB[varOffset + j];
          }
          component->y().setDataUpdated();
          component->yp().setDataUpdated();
          component->yB().setDataUpdated();
          component->ypB().setDataUpdated();
        }
        varOffset += component->size();

        if (component->sizeParams() > 0)
        {
          auto* component_param = component->param().getData();
          for (IdxT j = 0; j < component->sizeParams(); ++j)
          {
            component_param[j] = param[optOffset + j];
          }
          component->param().setDataUpdated();
        }
        optOffset += component->sizeParams();
      }

      for (const auto& bus : buses_)
      {
        bus->evaluateAdjointResidual();
      }

      for (const auto& component : components_)
      {
        component->evaluateAdjointResidual();
      }

      // Update residual vector
      IdxT resOffset = 0;
      for (const auto& bus : buses_)
      {
        if (bus->size() > 0)
        {
          const auto* bus_fB = bus->getAdjointResidual().getData();
          for (IdxT j = 0; j < bus->size(); ++j)
          {
            fB[resOffset + j] = bus_fB[j];
          }
        }
        resOffset += bus->size();
      }

      for (const auto& component : components_)
      {
        if (component->size() > 0)
        {
          const auto* component_fB = component->getAdjointResidual().getData();
          for (IdxT j = 0; j < component->size(); ++j)
          {
            fB[resOffset + j] = component_fB[j];
          }
        }
        resOffset += component->size();
      }

      fB_.setDataUpdated();

      return 0;
    }

    // int evaluateAdjointJacobian(){return 0;}

    /**
     * @brief Evaluate adjoint integrand for the system model.
     *
     * @pre Assumes there are no integrands in bus models.
     * @pre Assumes integrand is implemented in only _one_ component.
     *
     */
    int evaluateAdjointIntegrand()
    {
      // First, update variables
      IdxT           varOffset = 0;
      IdxT           optOffset = 0;
      const auto*    y         = y_.getData();
      const auto*    yp        = yp_.getData();
      const auto*    yB        = yB_.getData();
      const auto*    ypB       = ypB_.getData();
      const ScalarT* param     = nullptr;
      ScalarT*       gB        = nullptr;

      if (this->sizeParams() > 0)
      {
        param = param_.getData();
      }
      if (this->sizeQuadrature() > 0 && this->sizeParams() > 0)
      {
        gB = gB_.getData();
      }

      for (const auto& bus : buses_)
      {
        if (bus->size() > 0)
        {
          auto* bus_y   = bus->y().getData();
          auto* bus_yp  = bus->yp().getData();
          auto* bus_yB  = bus->yB().getData();
          auto* bus_ypB = bus->ypB().getData();

          for (IdxT j = 0; j < bus->size(); ++j)
          {
            bus_y[j]   = y[varOffset + j];
            bus_yp[j]  = yp[varOffset + j];
            bus_yB[j]  = yB[varOffset + j];
            bus_ypB[j] = ypB[varOffset + j];
          }
          bus->y().setDataUpdated();
          bus->yp().setDataUpdated();
          bus->yB().setDataUpdated();
          bus->ypB().setDataUpdated();
        }
        varOffset += bus->size();

        if (bus->sizeParams() > 0)
        {
          auto* bus_param = bus->param().getData();
          for (IdxT j = 0; j < bus->sizeParams(); ++j)
          {
            bus_param[j] = param[optOffset + j];
          }
          bus->param().setDataUpdated();
        }
        optOffset += bus->sizeParams();
      }

      for (const auto& component : components_)
      {
        if (component->size() > 0)
        {
          auto* component_y   = component->y().getData();
          auto* component_yp  = component->yp().getData();
          auto* component_yB  = component->yB().getData();
          auto* component_ypB = component->ypB().getData();

          for (IdxT j = 0; j < component->size(); ++j)
          {
            component_y[j]   = y[varOffset + j];
            component_yp[j]  = yp[varOffset + j];
            component_yB[j]  = yB[varOffset + j];
            component_ypB[j] = ypB[varOffset + j];
          }
          component->y().setDataUpdated();
          component->yp().setDataUpdated();
          component->yB().setDataUpdated();
          component->ypB().setDataUpdated();
        }
        varOffset += component->size();

        if (component->sizeParams() > 0)
        {
          auto* component_param = component->param().getData();
          for (IdxT j = 0; j < component->sizeParams(); ++j)
          {
            component_param[j] = param[optOffset + j];
          }
          component->param().setDataUpdated();
        }
        optOffset += component->sizeParams();
      }

      // Evaluate integrand and update global vector
      for (const auto& component : components_)
      {
        if (component->sizeQuadrature() == 1)
        {
          component->evaluateAdjointIntegrand();
          if (component->sizeParams() > 0)
          {
            const auto* component_gB = component->getAdjointIntegrand().getData();
            for (IdxT j = 0; j < size_opt_; ++j)
            {
              gB[j] = component_gB[j];
            }
          }
          break;
        }
      }
      if (this->sizeQuadrature() > 0 && this->sizeParams() > 0)
      {
        gB_.setDataUpdated();
      }
      return 0;
    }

    void updateTime(RealT t, RealT a)
    {
      for (const auto& component : components_)
      {
        component->updateTime(t, a);
      }
    }

    void addBus(bus_type* bus)
    {
      buses_.push_back(bus);
    }

    void addComponent(component_type* component)
    {
      components_.push_back(component);
    }

  private:
    std::vector<bus_type*>       buses_;
    std::vector<component_type*> components_;

  }; // class SystemModel

} // namespace GridKit
