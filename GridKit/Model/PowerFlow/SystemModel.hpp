
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

      auto allocate_host_vector = [](VectorT& vector, IdxT n)
      {
        vector.resize(n);
        vector.allocate(memory::HOST);
        vector.setDataUpdated(memory::HOST);
      };

      allocate_host_vector(yB_, size_);
      allocate_host_vector(ypB_, size_);
      allocate_host_vector(fB_, size_);
      allocate_host_vector(g_, size_quad_);
      allocate_host_vector(gB_, size_quad_ * size_opt_);
      allocate_host_vector(param_, size_opt_);
      allocate_host_vector(param_lo_, size_opt_);
      allocate_host_vector(param_up_, size_opt_);

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
      IdxT varOffset = 0;
      IdxT optOffset = 0;

      for (const auto& bus : buses_)
      {
        bus->initialize();
      }

      for (const auto& bus : buses_)
      {
        for (IdxT j = 0; j < bus->size(); ++j)
        {
          y_[varOffset + j]  = bus->y()[j];
          yp_[varOffset + j] = bus->yp()[j];
        }
        varOffset += bus->size();

        for (IdxT j = 0; j < bus->sizeParams(); ++j)
        {
          param_[optOffset + j]    = bus->param()[j];
          param_lo_[optOffset + j] = bus->param_lo()[j];
          param_up_[optOffset + j] = bus->param_up()[j];
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
        for (IdxT j = 0; j < component->size(); ++j)
        {
          y_[varOffset + j]  = component->y()[j];
          yp_[varOffset + j] = component->yp()[j];
        }
        varOffset += component->size();

        for (IdxT j = 0; j < component->sizeParams(); ++j)
        {
          param_[optOffset + j]    = component->param()[j];
          param_lo_[optOffset + j] = component->param_lo()[j];
          param_up_[optOffset + j] = component->param_up()[j];
        }
        optOffset += component->sizeParams();
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
      IdxT offset = 0;
      for (const auto& bus : buses_)
      {
        bus->setAbsoluteTolerance(rel_tol);
        for (IdxT j = 0; j < bus->size(); ++j)
        {
          abs_tol_[offset + j] = bus->absoluteTolerance()[j];
        }
        offset += bus->size();
      }

      for (const auto& component : components_)
      {
        component->setAbsoluteTolerance(rel_tol);
        for (IdxT j = 0; j < component->size(); ++j)
        {
          abs_tol_[offset + j] = component->absoluteTolerance()[j];
        }
        offset += component->size();
      }

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
      IdxT varOffset = 0;
      IdxT optOffset = 0;
      for (const auto& bus : buses_)
      {
        for (IdxT j = 0; j < bus->size(); ++j)
        {
          bus->y()[j]  = y_[varOffset + j];
          bus->yp()[j] = yp_[varOffset + j];
        }
        varOffset += bus->size();

        for (IdxT j = 0; j < bus->sizeParams(); ++j)
        {
          bus->param()[j] = param_[optOffset + j];
        }
        optOffset += bus->sizeParams();

        bus->evaluateResidual();
      }

      for (const auto& component : components_)
      {
        for (IdxT j = 0; j < component->size(); ++j)
        {
          component->y()[j]  = y_[varOffset + j];
          component->yp()[j] = yp_[varOffset + j];
        }
        varOffset += component->size();

        for (IdxT j = 0; j < component->sizeParams(); ++j)
        {
          component->param()[j] = param_[optOffset + j];
        }
        optOffset += component->sizeParams();

        component->evaluateResidual();
      }

      // Update residual vector
      IdxT resOffset = 0;
      for (const auto& bus : buses_)
      {
        for (IdxT j = 0; j < bus->size(); ++j)
        {
          f_[resOffset + j] = bus->getResidual()[j];
        }
        resOffset += bus->size();
      }

      for (const auto& component : components_)
      {
        for (IdxT j = 0; j < component->size(); ++j)
        {
          f_[resOffset + j] = component->getResidual()[j];
        }
        resOffset += component->size();
      }

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
      IdxT varOffset = 0;
      IdxT optOffset = 0;
      for (const auto& bus : buses_)
      {
        for (IdxT j = 0; j < bus->size(); ++j)
        {
          bus->y()[j]  = y_[varOffset + j];
          bus->yp()[j] = yp_[varOffset + j];
        }
        varOffset += bus->size();

        for (IdxT j = 0; j < bus->sizeParams(); ++j)
        {
          bus->param()[j] = param_[optOffset + j];
        }
        optOffset += bus->sizeParams();

        bus->evaluateIntegrand();
      }

      for (const auto& component : components_)
      {
        for (IdxT j = 0; j < component->size(); ++j)
        {
          component->y()[j]  = y_[varOffset + j];
          component->yp()[j] = yp_[varOffset + j];
        }
        varOffset += component->size();

        for (IdxT j = 0; j < component->sizeParams(); ++j)
        {
          component->param()[j] = param_[optOffset + j];
        }
        optOffset += component->sizeParams();

        component->evaluateIntegrand();
      }

      // Update integrand vector
      IdxT intOffset = 0;
      for (const auto& bus : buses_)
      {
        for (IdxT j = 0; j < bus->sizeQuadrature(); ++j)
        {
          g_[intOffset + j] = bus->getIntegrand()[j];
        }
        intOffset += bus->sizeQuadrature();
      }

      for (const auto& component : components_)
      {
        for (IdxT j = 0; j < component->sizeQuadrature(); ++j)
        {
          g_[intOffset + j] = component->getIntegrand()[j];
        }
        intOffset += component->sizeQuadrature();
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
      IdxT offset    = 0;
      IdxT optOffset = 0;

      // Update bus variables and optimization parameters
      for (const auto& bus : buses_)
      {
        for (IdxT j = 0; j < bus->size(); ++j)
        {
          bus->y()[j]  = y_[offset + j];
          bus->yp()[j] = yp_[offset + j];
        }
        offset += bus->size();

        for (IdxT j = 0; j < bus->sizeParams(); ++j)
        {
          bus->param()[j] = param_[optOffset + j];
        }
        optOffset += bus->sizeParams();
      }

      // Update component variables and optimization parameters
      for (const auto& component : components_)
      {
        for (IdxT j = 0; j < component->size(); ++j)
        {
          component->y()[j]  = y_[offset + j];
          component->yp()[j] = yp_[offset + j];
        }
        offset += component->size();

        for (IdxT j = 0; j < component->sizeParams(); ++j)
        {
          component->param()[j] = param_[optOffset + j];
        }
        optOffset += component->sizeParams();
      }

      // Reset counter
      offset = 0;

      // Initialize bus adjoints
      for (const auto& bus : buses_)
      {
        bus->initializeAdjoint();

        for (IdxT j = 0; j < bus->size(); ++j)
        {
          yB_[offset + j]  = bus->yB()[j];
          ypB_[offset + j] = bus->ypB()[j];
        }
        offset += bus->size();
      }

      // Initialize component adjoints
      for (const auto& component : components_)
      {
        component->initializeAdjoint();

        for (IdxT j = 0; j < component->size(); ++j)
        {
          yB_[offset + j]  = component->yB()[j];
          ypB_[offset + j] = component->ypB()[j];
        }
        offset += component->size();
      }

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
      IdxT varOffset = 0;
      IdxT optOffset = 0;

      // Update variables in component models
      for (const auto& bus : buses_)
      {
        for (IdxT j = 0; j < bus->size(); ++j)
        {
          bus->y()[j]   = y_[varOffset + j];
          bus->yp()[j]  = yp_[varOffset + j];
          bus->yB()[j]  = yB_[varOffset + j];
          bus->ypB()[j] = ypB_[varOffset + j];
        }
        varOffset += bus->size();

        for (IdxT j = 0; j < bus->sizeParams(); ++j)
        {
          bus->param()[j] = param_[optOffset + j];
        }
        optOffset += bus->sizeParams();
      }

      for (const auto& component : components_)
      {
        for (IdxT j = 0; j < component->size(); ++j)
        {
          component->y()[j]   = y_[varOffset + j];
          component->yp()[j]  = yp_[varOffset + j];
          component->yB()[j]  = yB_[varOffset + j];
          component->ypB()[j] = ypB_[varOffset + j];
        }
        varOffset += component->size();

        for (IdxT j = 0; j < component->sizeParams(); ++j)
        {
          component->param()[j] = param_[optOffset + j];
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
        for (IdxT j = 0; j < bus->size(); ++j)
        {
          fB_[resOffset + j] = bus->getAdjointResidual()[j];
        }
        resOffset += bus->size();
      }

      for (const auto& component : components_)
      {
        for (IdxT j = 0; j < component->size(); ++j)
        {
          fB_[resOffset + j] = component->getAdjointResidual()[j];
        }
        resOffset += component->size();
      }

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
      IdxT varOffset = 0;
      IdxT optOffset = 0;
      for (const auto& bus : buses_)
      {
        for (IdxT j = 0; j < bus->size(); ++j)
        {
          bus->y()[j]   = y_[varOffset + j];
          bus->yp()[j]  = yp_[varOffset + j];
          bus->yB()[j]  = yB_[varOffset + j];
          bus->ypB()[j] = ypB_[varOffset + j];
        }
        varOffset += bus->size();

        for (IdxT j = 0; j < bus->sizeParams(); ++j)
        {
          bus->param()[j] = param_[optOffset + j];
        }
        optOffset += bus->sizeParams();
      }

      for (const auto& component : components_)
      {
        for (IdxT j = 0; j < component->size(); ++j)
        {
          component->y()[j]   = y_[varOffset + j];
          component->yp()[j]  = yp_[varOffset + j];
          component->yB()[j]  = yB_[varOffset + j];
          component->ypB()[j] = ypB_[varOffset + j];
        }
        varOffset += component->size();

        for (IdxT j = 0; j < component->sizeParams(); ++j)
        {
          component->param()[j] = param_[optOffset + j];
        }
        optOffset += component->sizeParams();
      }

      // Evaluate integrand and update global vector
      for (const auto& component : components_)
      {
        if (component->sizeQuadrature() == 1)
        {
          component->evaluateAdjointIntegrand();
          for (IdxT j = 0; j < size_opt_; ++j)
          {
            gB_[j] = component->getAdjointIntegrand()[j];
          }
          break;
        }
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
