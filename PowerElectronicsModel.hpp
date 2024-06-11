

#pragma once

#include <iostream>
#include <vector>
#include <cassert>

#include <ScalarTraits.hpp>
#include <ModelEvaluatorImpl.hpp>
#include <CircuitGraph.hpp>
#include <ComponentLib/PowerElectronics/CircuitComponent.hpp>

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
         * @brief Default constructor for the system model
         */
        PowerElectronicsModel() : ModelEvaluatorImpl<ScalarT, IdxT>(0, 0, 0)
        {
            // Set system model tolerances as default
            rtol_ = 1e-4;
            atol_ = 1e-4;
            this->max_steps_ = 2000;
            // By default not use jacobian
            usejac_ = false;
        }

        /**
         * @brief Constructor for the system model
         */
        PowerElectronicsModel(double rt = 1e-4, double at = 1e-4, bool ju = false, IdxT msa = 2000) : ModelEvaluatorImpl<ScalarT, IdxT>(0, 0, 0)
        {
            // Set system model tolerances from input
            rtol_ = rt;
            atol_ = at;
            this->max_steps_ = msa;
            // Can choose if to use jacobain
            usejac_ = ju;
        }

        /**
         * @brief Destructor for the system model
         */
        virtual ~PowerElectronicsModel()
        {
            for (auto comp : this->components_)
                delete comp;
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
         * @brief Will check if each component has jacobian avalible. If one doesn't then jacobain is false.
         *
         *
         * @return true
         * @return false
         */
        bool hasJacobian()
        {
            if (!this->usejac_)
                return false;

            for (const auto &component : components_)
            {
                if (!component->hasJacobian())
                {
                    return false;
                }
            }
            return true;
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
            size_ = s;
            for (const auto &component : components_)
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
            for (const auto &component : components_)
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
            for (const auto &component : components_)
            {
                for (IdxT j = 0; j < component->size(); ++j)
                {
                    if (component->getNodeConnection(j) != neg1_)
                    {
                        component->y()[j] = y_[component->getNodeConnection(j)];
                        component->yp()[j] = yp_[component->getNodeConnection(j)];
                    }
                    else
                    {
                        component->y()[j] = 0.0;
                        component->yp()[j] = 0.0;
                    }
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
                f_[i] = 0.0;
            }

            this->distributeVectors();

            // Update system residual vector

            for (const auto &component : components_)
            {
                // TODO:check return type
                component->evaluateResidual();
                for (IdxT j = 0; j < component->size(); ++j)
                {
                    //@todo should do a different grounding check
                    if (component->getNodeConnection(j) != neg1_)
                    {
                        f_[component->getNodeConnection(j)] += component->getResidual()[j];
                    }
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
            J_.zeroMatrix();
            distributeVectors();

            // Evaluate component jacs
            for (const auto &component : components_)
            {
                component->evaluateJacobian();

                // get references to local jacobain
                std::tuple<std::vector<IdxT>&, std::vector<IdxT>&, std::vector<ScalarT>&> tpm = component->getJacobian().getEntries();
                const auto& [r, c, v] = tpm;

                // Create copies of data to handle groundings
                std::vector<IdxT> rgr;
                std::vector<IdxT> cgr;
                std::vector<ScalarT> vgr;
                for (IdxT i = 0; i < static_cast<IdxT>(r.size()); i++)
                {
                    if (component->getNodeConnection(r[i]) != neg1_ && component->getNodeConnection(c[i]) != neg1_)
                    {
                        rgr.push_back(component->getNodeConnection(r[i]));
                        cgr.push_back(component->getNodeConnection(c[i]));
                        vgr.push_back(v[i]);
                    }
                }

                // AXPY to Global Jacobian
                J_.axpy(1.0, rgr, cgr, vgr);
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
            for (const auto &component : components_)
            {
                component->updateTime(t, a);
            }
            time_ = t;
            alpha_ = a;
        }

        void addComponent(component_type *component)
        {
            components_.push_back(component);
        }

    private:

        static constexpr IdxT neg1_ = static_cast<IdxT>(-1);

        std::vector<component_type*> components_;
        bool usejac_;

    }; // class PowerElectronicsModel

} // namespace ModelLib
