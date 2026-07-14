
/**
 * @file SystemSteadyStaeModel.hpp
 * @author Slaven Peles <slaven.peles@pnnl.gov>
 *
 * Contains definition of power flow analysis class.
 *
 */
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
   * ModelEvaluator for the system model. This is still work in
   * progress and code is not optimized.
   *
   * @todo Address thread safety for the system model methods.
   *
   * @todo Tolerance management needs to be reconsidered.
   *
   */
  template <class ScalarT, typename IdxT>
  class SystemSteadyStateModel : public ModelEvaluatorImpl<ScalarT, IdxT>
  {
    using bus_type       = BaseBus<ScalarT, IdxT>;
    using component_type = ModelEvaluatorImpl<ScalarT, IdxT>;
    using RealT          = typename ModelEvaluatorImpl<ScalarT, IdxT>::RealT;

    using ModelEvaluatorImpl<ScalarT, IdxT>::size_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::size_quad_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::size_opt_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::nnz_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::time_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::alpha_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::y_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::yp_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::yB_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::ypB_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::tag_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::abs_tol_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::f_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::fB_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::g_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::gB_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::param_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::param_up_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::param_lo_;

  public:
    /**
     * @brief Constructor for the system model
     */
    SystemSteadyStateModel()
      : ModelEvaluatorImpl<ScalarT, IdxT>(0, 0, 0)
    {
    }

    /**
     * @brief Construct a new System Steady State Model object. Allows for simple allocation.
     *
     * @param mp model data
     */
    SystemSteadyStateModel(GridKit::PowerFlowData::SystemModelData<ScalarT, IdxT> mp)
      : ModelEvaluatorImpl<ScalarT, IdxT>(0, 0, 0)
    {
      // add buses
      for (auto busdata : mp.bus)
      {
        auto* bus = BusFactory<double, size_t>::create(busdata);
        this->addBus(bus);
      }

      // add generators
      for (auto gendata : mp.gen)
      {
        auto* gen = GeneratorFactory<double, size_t>::create(this->getBus(gendata.bus), gendata);
        this->addComponent(gen);
      }

      // add branches
      for (auto branchdata : mp.branch)
      {
        auto* branch = new Branch<double, size_t>(this->getBus(branchdata.fbus), this->getBus(branchdata.tbus), branchdata);
        this->addComponent(branch);
      }

      // add loads
      for (auto loaddata : mp.load)
      {
        auto* loadm = new Load<double, size_t>(this->getBus(loaddata.bus_i), loaddata);
        this->addComponent(loadm);
      }

      // There is no Generator Cost Object
      // TODO: Implment for GenCost
    }

    /**
     * @brief Destructor for the system model
     */
    virtual ~SystemSteadyStateModel()
    {
      for (auto comp : this->components_)
        delete comp;
      for (auto bus : this->buses_)
        delete bus;
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
      size_ = 0;

      // Allocate all buses
      for (const auto& bus : buses_)
      {
        bus->allocate();
        size_ += bus->size();
      }

      // Allocate all components
      for (const auto& component : components_)
      {
        component->allocate();
        size_ += component->size();
      }

      // Allocate global vectors
      this->allocateVectors(size_);
      tag_.resize(size_);

      return 0;
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
      IdxT  varOffset = 0;
      auto* y         = y_.getData();

      for (const auto& bus : buses_)
      {
        bus->initialize();
      }

      for (const auto& bus : buses_)
      {
        if (bus->size() > 0)
        {
          const auto* bus_y = bus->y().getData();

          for (IdxT j = 0; j < bus->size(); ++j)
          {
            y[varOffset + j] = bus_y[j];
          }
        }
        varOffset += bus->size();
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
          const auto* component_y = component->y().getData();

          for (IdxT j = 0; j < component->size(); ++j)
          {
            y[varOffset + j] = component_y[j];
          }
        }
        varOffset += component->size();
      }
      y_.setDataUpdated();
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
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
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
      IdxT        varOffset = 0;
      const auto* y         = y_.getData();
      auto*       f         = f_.getData();

      for (const auto& bus : buses_)
      {
        if (bus->size() > 0)
        {
          auto* bus_y = bus->y().getData();

          for (IdxT j = 0; j < bus->size(); ++j)
          {
            bus_y[j] = y[varOffset + j];
          }
          bus->y().setDataUpdated();
        }
        varOffset += bus->size();
        bus->evaluateResidual();
      }

      for (const auto& component : components_)
      {
        if (component->size() > 0)
        {
          auto* component_y = component->y().getData();

          for (IdxT j = 0; j < component->size(); ++j)
          {
            component_y[j] = y[varOffset + j];
          }
          component->y().setDataUpdated();
        }
        varOffset += component->size();
        component->evaluateResidual();
      }

      // Update system residual vector
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
      return 0;
    }

    void updateTime(RealT /* t */, RealT /* a */)
    {
    }

    void addBus(bus_type* bus)
    {
      buses_.push_back(bus);
    }

    void addComponent(component_type* component)
    {
      components_.push_back(component);
    }

    bus_type* getBus(IdxT busid)
    {
      // Need to implement mapping of bus IDs to buses in the system model
      assert((buses_[busid - 1])->BusID() == busid);
      return buses_[busid - 1];
    }

  private:
    std::vector<bus_type*>       buses_;
    std::vector<component_type*> components_;

  }; // class SystemSteadyStateModel

} // namespace GridKit
