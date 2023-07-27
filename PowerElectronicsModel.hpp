

#pragma once

#include <iostream>
#include <vector>
#include <cassert>

#include <ScalarTraits.hpp>
#include <ModelEvaluatorImpl.hpp>
#include <CircuitGraph.hpp>
#include <ComponentLib/PowerElectronicsComponents/CircuitComponent.hpp>

namespace ModelLib
{

template <class ScalarT, typename IdxT>
class PowerElectronicsModel : public ModelEvaluatorImpl<ScalarT, IdxT>
{
    typedef CircuitComponent<ScalarT, IdxT> component_type;

    using ModelEvaluatorImpl<ScalarT, IdxT>::size_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::size_quad_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::size_opt_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::nnz_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::time_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::alpha_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::y_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::yp_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::yB_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::ypB_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::tag_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::f_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::fB_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::g_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::gB_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::rtol_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::atol_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::param_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::param_up_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::param_lo_;

public:
    /**
     * @brief Constructor for the system model
     */
    PowerElectronicsModel() : ModelEvaluatorImpl<ScalarT, IdxT>(0, 0, 0)
    {
        // Set system model tolerances
        rtol_ = 1e-5;
        atol_ = 1e-5;
    }

    /**
     * @brief Destructor for the system model
     */
    virtual ~PowerElectronicsModel()
    {
        for (auto comp : this->components_) delete comp;
    }

    /**
     * @brief allocator default
     * 
     * @todo this should throw an exception as no allocation without a graph is allowed. Or needs to be removed from the base class
     * 
     * @return int 
     */
    int allocate()
    {
        
        return 1;
    }

    /**
     * @brief Allocate the vector data with size amount
     * @todo Add capability to go through component model connection to get the size of the actual vector
     * 
     * @param s 
     * @return int 
     */
    int allocate(IdxT s)
    {

        // Allocate all components
        this->size_ = s;
        for(const auto& component : components_)
        {
            component->allocate();
        }

        // Allocate global vectors
        y_.resize(size_);
        yp_.resize(size_);
        f_.resize(size_);

        return 0;
    }

    int initialize()
    {

        // Initialize components
        for(const auto& component : components_)
        {
            component->initialize();
        }

        for(const auto& component : components_)
        {
            for(IdxT j=0; j<component->size(); ++j)
            {
                component->y()[j] = y_[component->getNodeConnection(j)];
                component->yp()[j] = yp_[component->getNodeConnection(j)];
            }
        }
        return 0;
    }

    int tagDifferentiable()
    {
        return 0;
    }

    int evaluateResidual()
    {
        for (IdxT i = 0; i < this->f_.size(); i++)
        {
            f_[i] = 0;
        }
        
        // Update variables
        for(const auto& component : components_)
        {
            for(IdxT j=0; j<component->size(); ++j)
            {
                component->y()[j]  = y_[component->getNodeConnection(j)];
				component->yp()[j] = yp_[component->getNodeConnection(j)];
            }
            component->evaluateResidual();
        }

        // Update system residual vector

        for(const auto& component : components_)
        {
            for(IdxT j=0; j<component->size(); ++j)
            {
                f_[component->getNodeConnection(j)] += component->getResidual()[j];
            }
        }

        return 0;
    }

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
     *
     */
    int evaluateAdjointResidual()
    {
        return 0;
    }


    /**
     * @brief Evaluate adjoint integrand for the system model.
     *
     *
     */
    int evaluateAdjointIntegrand()
    {
        return 0;
    }

    void updateTime(ScalarT t, ScalarT a)
    {
    }

    void addComponent(component_type* component)
    {
        this->components_.push_back(component);
    }

private:
    std::vector<component_type*> components_;

}; // class PowerElectronicsModel

} // namespace ModelLib
