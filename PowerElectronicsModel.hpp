

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
    using ModelEvaluatorImpl<ScalarT, IdxT>::time_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::alpha_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::y_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::yp_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::yB_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::ypB_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::tag_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::f_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::fB_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::g_;
    // using ModelEvaluatorImpl<ScalarT, IdxT>::gB_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::J_;
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

    /**
     * @brief Set intial y and y' of each component
     * 
     * @return int 
     */
    int initialize()
    {

        // Initialize components
        for(const auto& component : components_)
        {
            component->initialize();
        }
        this->distributeVectors();
        
        return 0;
    }

    /**
     * @brief Distribute y and y' to each component based of node connection graph
     * 
     * @return int 
     */
    int distributeVectors()
    {
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

    /**
     * @brief Evaluate Residuals at each component then collect them
     * 
     * @return int 
     */
    int evaluateResidual()
    {
        for (IdxT i = 0; i < this->f_.size(); i++)
        {
            f_[i] = 0;
        }
        
        this->distributeVectors();

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

    /**
     * @brief Creates the Sparse COO Jacobian representing  \alpha dF/dy' + dF/dy
     * 
     * @return int 
     */
    int evaluateJacobian()
	{
        this->J_.zeroMatrix();
        this->distributeVectors();

        //Evaluate component jacs
        for(const auto& component : components_)
        {
            component->evaluateJacobian();
            std::vector<IdxT> r;
            std::vector<IdxT> c;
            std::vector<ScalarT> v;
	        std::tie(r, c, v) = component->getJacobian().getEntrieCopies();
            for (IdxT i = 0; i < static_cast<IdxT>(r.size()); i++)
            {
                r[i] = component->getNodeConnection(r[i]);
                c[i] = component->getNodeConnection(c[i]);
            }
            this->J_.AXPY(1.0, r, c, v);
        }
        
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

    /**
     * @brief Distribute time and time scaling for each component
     * 
     * @param t 
     * @param a 
     */
    void updateTime(ScalarT t, ScalarT a)
    {
        for(const auto& component : components_)
        {
            component->updateTime(t, a);
        }
    }

    void addComponent(component_type* component)
    {
        this->components_.push_back(component);
    }

private:
    std::vector<component_type*> components_;

}; // class PowerElectronicsModel

} // namespace ModelLib
