

#pragma once

#include <cassert>
#include <unordered_map>
#include <vector>

#include <GridKit/Constants.hpp>
#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>
#include <GridKit/Model/PowerElectronics/CircuitNode.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridBusDQ/MicrogridBusDQ.hpp>
#include <GridKit/ScalarTraits.hpp>

#include "GridKit/Model/PowerElectronics/MicrogridBusDQ/MicrogridBusDQ.hpp"

namespace GridKit
{
  template <class ScalarT, typename IdxT>
  class SubsystemModel : public CircuitComponent<ScalarT, IdxT>
  {
    using RealT          = typename CircuitComponent<ScalarT, IdxT>::RealT;
    using MatrixT        = typename CircuitComponent<ScalarT, IdxT>::MatrixT;
    using CsrMatrixT     = typename CircuitComponent<ScalarT, IdxT>::CsrMatrixT;
    using component_type = CircuitComponent<ScalarT, IdxT>;
    using node_type      = CircuitNode<ScalarT, IdxT>;

    using CircuitComponent<ScalarT, IdxT>::size_;
    using CircuitComponent<ScalarT, IdxT>::n_extern_;
    using CircuitComponent<ScalarT, IdxT>::n_intern_;
    using CircuitComponent<ScalarT, IdxT>::nnz_;
    using CircuitComponent<ScalarT, IdxT>::time_;
    using CircuitComponent<ScalarT, IdxT>::alpha_;
    using CircuitComponent<ScalarT, IdxT>::y_;
    using CircuitComponent<ScalarT, IdxT>::yp_;
    using CircuitComponent<ScalarT, IdxT>::f_;
    using CircuitComponent<ScalarT, IdxT>::jac_;
    using CircuitComponent<ScalarT, IdxT>::rel_tol_;
    using CircuitComponent<ScalarT, IdxT>::abs_tol_;

  public:
    /**
     * @brief Default constructor for the system model
     *
     * @post System model parameters set as default
     */
    SubsystemModel(IdxT refframe_index = neg1_)
    {
      // Set system model parameters as default
      rel_tol_         = 1e-4;
      abs_tol_         = 1e-4;
      this->max_steps_ = 2000;
      // By default don't use the jacobian
      use_jac_         = false;
      refframe_index_  = refframe_index;
    }

    virtual ~SubsystemModel()
    {
      for (auto comp : this->components_)
      {
        delete comp;
      }
      delete csr_jac_;
      delete[] map_to_csr_;
    }

    bool hasJacobian() final
    {
      if (!this->use_jac_)
        return false;

      for (const auto& component : components_)
      {
        if (!component->hasJacobian())
        {
          return false;
        }
      }
      return true;
    }

    int allocate() override
    {
      size_t counter = 0;

      // First pass: Map global component indices to local subsystem indices.
      for (const auto& component : components_)
      {
        component->allocate();

        for (IdxT i = component->getExternSize(); i < component->size(); i++)
        {

          IdxT index = component->getNodeConnection(i);

          if (index != neg1_)
          {
            internal_map_[index] = counter;
            counter++;
          }
        }

        // TODO : To be removed once MicrogridBusDQ owns its internals
        auto* comp = dynamic_cast<GridKit::MicrogridBusDQ<double, size_t>*>(component);

        if (comp != nullptr)
        {
          IdxT index1 = component->getNodeConnection(0);
          IdxT index2 = component->getNodeConnection(1);

          internal_map_[index1] = counter;
          counter++;
          internal_map_[index2] = counter;
          counter++;
        }
      }

      // If this subsystem contains the reference generator, add the rotor motor as an internal
      if (refframe_index_ != neg1_)
      {
        internal_map_[refframe_index_] = counter;
        counter++;
      }

      // Second pass: Identify variables required by this partition but owned by another partition.
      for (const auto& component : components_)
      {
        for (IdxT i = 0; i < component->getExternSize(); i++)
        {
          IdxT index = component->getNodeConnection(i);

          if (internal_map_.count(index) < 1 && external_map_.count(index) < 1 && index != neg1_)
          {
            external_map_[index] = counter;
            counter++;
          }
        }
      }

      n_intern_ = internal_map_.size();
      n_extern_ = external_map_.size();
      size_     = n_intern_ + n_extern_;

      this->mapGlobalToLocal();

      // Allocate subsystem vectors.
      y_.resize(n_intern_);
      yp_.resize(n_intern_);
      f_.resize(n_intern_);

      // Allocate vectors associated with external coupling variables.
      y_ext_.resize(n_extern_);
      yp_ext_.resize(n_extern_);
      f_ext_.resize(n_extern_);
      external_indices_.resize(n_extern_);

      // Store the mapping from local subsystem indices back to their global system indices
      for (const auto& [global_idx, local_idx] : internal_map_)
      {
        this->setExternalConnectionNodes(local_idx, global_idx);
      }

      // Store the global indices of all external coupling variables.
      for (const auto& [global_idx, local_idx] : external_map_)
      {
        external_indices_[local_idx - n_intern_] = global_idx;
      }

      return 0;
    }

    /**
     * @brief Set intial y and y' of each component
     *
     * @return int 0 if successful, positive if there's a recoverable error, negative if unrecoverable
     */
    int initialize() final
    {
      // Initialize components
      for (const auto& component : components_)
      {
        component->initialize();
      }
      this->distributeVectors();

      return 0;
    }

    /**
     * @brief Distribute y and y' to each component based of node connection graph
     *
     * @post Each component has y and y' set
     *
     * @return int 0 if successful, positive if there's a recoverable error, negative if unrecoverable
     */
    int distributeVectors()
    {
      for (const auto& component : components_)
      {
        IdxT                  size = component->size();
        std::vector<ScalarT>& y    = component->y();
        std::vector<ScalarT>& yp   = component->yp();

        for (IdxT j = 0; j < size; ++j)
        {
          if (component->getNodeConnection(j) == neg1_)
          {
            y[j]  = 0.0;
            yp[j] = 0.0;
          }
          else if (component->getNodeConnection(j) < this->getInternalSize())
          {
            y[j]  = y_[component->getNodeConnection(j)];
            yp[j] = yp_[component->getNodeConnection(j)];
          }
          else
          {
            y[j]  = y_ext_[component->getNodeConnection(j) - this->getInternalSize()];
            yp[j] = yp_ext_[component->getNodeConnection(j) - this->getInternalSize()];
          }
        }
      }
      return 0;
    }

    int tagDifferentiable() final
    {
      return 0;
    }

    /**
     * @brief Evaluate Residuals at each component then collect them
     *
     * @return int 0 if successful, positive if there's a recoverable error, negative if unrecoverable
     */
    int evaluateResidual() final
    {
      for (IdxT i = 0; i < this->f_.size(); i++)
      {
        f_[i] = 0.0;
      }

      this->distributeVectors();

      // Update system residual vector
      for (const auto& component : components_)
      {
        // TODO:check return type
        component->evaluateResidual();

        IdxT                        size     = component->size();
        const std::vector<ScalarT>& residual = component->getResidual();
        for (IdxT j = 0; j < size; ++j)
        {
          //@todo should do a different grounding check
          if (component->getNodeConnection(j) != neg1_ && component->getNodeConnection(j) < this->getInternalSize())
          {
            f_[component->getNodeConnection(j)] += residual[j];
          }
        }
      }

      return 0;
    }

    /**
     * @brief Creates the system Jacobian representing \alpha dF/dy' + dF/dy
     *
     *
     * @return int 0 if successful, positive if there's a recoverable error, negative if unrecoverable
     */
    int evaluateJacobian() final
    {
      return 0;
    }

    /**
     * @brief Evaluate integrands for the system quadratures.
     */
    int evaluateIntegrand() final
    {
      return 0;
    }

    /**
     * @brief Initialize system adjoint.
     *
     * Updates variables and optimization parameters, then initializes
     * adjoints locally and copies them to the system adjoint vector.
     */
    int initializeAdjoint() final
    {
      return 0;
    }

    /**
     * @brief Compute adjoint residual for the system model.
     *
     *
     */
    int evaluateAdjointResidual() final
    {
      return 0;
    }

    /**
     * @brief Evaluate adjoint integrand for the system model.
     *
     *
     */
    int evaluateAdjointIntegrand() final
    {
      return 0;
    }

    /**
     * @brief Distribute time and time scaling for each component
     *
     * @param t
     * @param a
     */
    void updateTime(RealT t, RealT a) final
    {
      for (const auto& component : components_)
      {
        component->updateTime(t, a);
      }
      time_  = t;
      alpha_ = a;
    }

    /**
     * @brief print the system residual in COO format
     *
     * @param[in] filename
     * @param[in] title
     */
    void printResidualMatrixMarket(std::string filename, std::string title)
    {
      writeVectorToMatrixMarket(f_, filename, title);
    }

    /**
     * @brief print the system Jacobian in COO format
     *
     * @param[in] filename
     * @param[in] title
     */
    void printJacobianMatrixMarket(std::string filename, std::string title)
    {
      jac_.printMatrixMarket(filename, title);
    }

    CsrMatrixT* getCsrJacobian() const override
    {
      return csr_jac_;
    }

    void addComponent(component_type* component)
    {
      components_.push_back(component);
    }

    std::unordered_map<IdxT, IdxT>& getInternalMap()
    {
      return internal_map_;
    }

    std::unordered_map<IdxT, IdxT>& getExternalMap()
    {
      return external_map_;
    }

    int mapGlobalToLocal()
    {
      for (const auto& component : components_)
      {

        for (IdxT i = 0; i < component->size(); i++)
        {
          IdxT index = component->getNodeConnection(i);

          if (index == neg1_)
          {
            continue;
          }
          else if (internal_map_.count(index) > 0)
          {
            component->setExternalConnectionNodes(i, internal_map_.at(index));
          }
          else
          {
            component->setExternalConnectionNodes(i, external_map_.at(index));
          }
        }
      }

      return 0;
    }

    int mapLocalToGlobal()
    {
      return 0;
    }

    std::vector<RealT>& getExternalDataY()
    {
      return y_ext_;
    }

    std::vector<RealT>& getExternalDataYP()
    {
      return yp_ext_;
    }

    std::vector<IdxT>& getExternalIndices()
    {
      return external_indices_;
    }

  private:
    static constexpr IdxT neg1_ = INVALID_INDEX<IdxT>;

    std::vector<component_type*>   components_;
    std::unordered_map<IdxT, IdxT> internal_map_;
    std::unordered_map<IdxT, IdxT> external_map_;
    std::vector<IdxT>              external_indices_;
    std::vector<RealT>             y_ext_;
    std::vector<RealT>             yp_ext_;
    std::vector<RealT>             f_ext_;
    IdxT                           refframe_index_{neg1_};

    IdxT*       map_to_csr_{nullptr};
    CsrMatrixT* csr_jac_{nullptr};

    int  jac_call_count_{0};
    bool use_jac_;

  }; // class SubsystemModel

} // namespace GridKit
