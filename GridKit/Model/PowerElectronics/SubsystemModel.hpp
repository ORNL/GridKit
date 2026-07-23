

#pragma once

#include <cassert>
#include <functional>
#include <optional>
#include <unordered_map>
#include <vector>

#include <GridKit/Constants.hpp>
#include <GridKit/Model/PowerElectronics/Bus/MicrogridBus.hpp>
#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>
#include <GridKit/Model/PowerElectronics/CircuitNode.hpp>
#include <GridKit/Model/PowerElectronics/SystemModelPowerElectronics.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{

  template <class ScalarT, typename IdxT>
  class SubsystemModel : public PowerElectronicsModel<ScalarT, IdxT>
  {

    using SystemModel    = PowerElectronicsModel<ScalarT, IdxT>;
    using RealT          = typename CircuitComponent<ScalarT, IdxT>::RealT;
    using CsrMatrixT     = typename CircuitComponent<ScalarT, IdxT>::CsrMatrixT;
    using component_type = CircuitComponent<ScalarT, IdxT>;
    using node_type      = PowerElectronics::NodeBase<ScalarT, IdxT>;
    using ForcingData    = std::pair<std::vector<ScalarT>, std::vector<ScalarT>>;
    using TimeFunction   = std::function<ForcingData(ScalarT)>;

    using SystemModel::abs_tol_;
    using SystemModel::allocated_;
    using SystemModel::allocateVectors;
    using SystemModel::alpha_;
    using SystemModel::f_;
    using SystemModel::f_int_;
    using SystemModel::n_extern_;
    using SystemModel::n_intern_;
    using SystemModel::nnz_;
    using SystemModel::size_;
    using SystemModel::tag_;
    using SystemModel::time_;
    using SystemModel::y_;
    using SystemModel::y_int_;
    using SystemModel::yp_;
    using SystemModel::yp_int_;

    using SystemModel::components_;
    using SystemModel::csr_jac_;
    using SystemModel::jac_call_count_;
    using SystemModel::map_to_csr_;
    using SystemModel::neg1_;
    using SystemModel::nodes_;
    using SystemModel::use_jac_;

  public:
    /**
     * @brief Default constructor for the system model
     *
     * @post System model parameters set as default
     */

    SubsystemModel(bool use_jac = false) : SystemModel(use_jac)
    {
    }

    ~SubsystemModel() override
    {
      // SubsystemModel does not own components_/nodes_ — they belong to (and are
      // deleted by) the parent system model.
      components_.clear();
      nodes_.clear();
    }

    int allocate() override
    {

      this->createGlobalToInternalMap();

      this->mapGlobalToLocal();

      n_intern_ = internal_map_.size();
      n_extern_ = external_map_.size();
      size_     = n_intern_;

      CircuitComponent<ScalarT, IdxT>::allocate();

      // Allocation always rebuilds the system Jacobian and its COO-to-CSR map.
      delete csr_jac_;
      csr_jac_ = nullptr;

      delete[] map_to_csr_;
      map_to_csr_ = nullptr;

      tag_.resize(size_);

      // Allocate subsystem vectors.
      y_ext_.resize(n_extern_);
      yp_ext_.resize(n_extern_);
      f_ext_.resize(n_extern_);

      external_indices_.resize(n_extern_);

      // Store the mapping from local subsystem indices back to their global system indices
      for (const auto [global_idx, local_idx] : internal_map_)
      {
        this->setExternalConnectionNodes(local_idx, global_idx);
      }

      // Store the global indices of all external coupling variables.
      for (const auto [global_idx, local_idx] : external_map_)
      {
        external_indices_[local_idx - n_intern_] = global_idx;
      }

      size_t      component_internal_idx = 0;
      const auto* y                      = y_.getData();
      const auto* yp                     = yp_.getData();
      auto*       f                      = f_.getData();
      for (component_type* comp : components_)
      {

        comp->setInternalPointer(&y[component_internal_idx]);
        comp->setInternalDerivativePointer(&yp[component_internal_idx]);
        comp->setInternalResidualPointer(&f[component_internal_idx]);

        component_internal_idx += comp->getInternalSize();
      }

      // Evaluate component Jacobians to get sparsity
      distributeVectors();
      for (component_type* component : components_)
      {
        component->evaluateJacobian();
      }

      auto isValidEntry = [this](IdxT row, IdxT col)
      {
        if (row == neg1_ || col == neg1_)
        {
          return false;
        }

        const bool row_is_internal = row < this->getInternalSize();
        const bool col_is_internal = col < this->getInternalSize();

        return (row_is_internal && col_is_internal);
      };

      IdxT nnz_dup = 0;

      for (const component_type* component : components_)
      {
        const IdxT* r   = component->jacobianCooRows();
        const IdxT* c   = component->jacobianCooCols();
        const IdxT  nnz = component->nnz();

        for (IdxT i = 0; i < nnz; ++i)
        {
          const IdxT row = component->getNodeConnection(r[i]);
          const IdxT col = component->getNodeConnection(c[i]);

          if (isValidEntry(row, col))
          {
            ++nnz_dup;
          }
        }
      }

      // Allocate COO triplet arrays (we own these until we hand off to CsrMatrix)
      IdxT*  rows_dup = new IdxT[nnz_dup];
      IdxT*  cols_dup = new IdxT[nnz_dup];
      RealT* vals_dup = new RealT[nnz_dup];

      IdxT counter = 0;

      for (const component_type* component : components_)
      {
        const IdxT*  r   = component->jacobianCooRows();
        const IdxT*  c   = component->jacobianCooCols();
        const RealT* v   = component->jacobianCooValues();
        const IdxT   nnz = component->nnz();

        for (IdxT i = 0; i < nnz; ++i)
        {
          const IdxT row = component->getNodeConnection(r[i]);
          const IdxT col = component->getNodeConnection(c[i]);

          if (!isValidEntry(row, col))
          {
            continue;
          }

          rows_dup[counter] = row;
          cols_dup[counter] = col;
          vals_dup[counter] = v[i];

          ++counter;
        }
      }

      // Build the system COO Jacobian
      LinearAlgebra::CooMatrix<RealT, IdxT> jac(size_, size_, nnz_dup, &rows_dup, &cols_dup, &vals_dup);

      // Populate CSR data with sort and deduplicate
      IdxT* row_ptrs = jac.getCsrRowData();

      // Deduplicated nnz
      nnz_ = jac.getNnz();

      // Allocate cols/vals with deduplicated nnz
      IdxT*  cols = new IdxT[nnz_];
      RealT* vals = new RealT[nnz_];

      std::copy(jac.getColData(), jac.getColData() + nnz_, cols);
      std::copy(jac.getValues(), jac.getValues() + nnz_, vals);

      // Create the CSR Jacobian
      csr_jac_ = new CsrMatrixT(size_, size_, nnz_, &row_ptrs, &cols, &vals);

      const IdxT* map_to_sorted = jac.getMapToSorted();
      const IdxT* map_to_dedup  = jac.getMapToDeduplicated();

      // Build a mappping from original COO index to CSR index
      map_to_csr_ = new IdxT[nnz_dup];
      for (IdxT i = 0; i < nnz_dup; ++i)
      {
        map_to_csr_[map_to_sorted[i]] = map_to_dedup[i];
      }

      allocated_ = true;
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

      if (forcing_function_)
      {
        const auto forcing = (*forcing_function_)(time_);

        const auto& y_forcing  = forcing.first;
        const auto& yp_forcing = forcing.second;

        assert(y_forcing.size() == y_ext_.size());
        assert(yp_forcing.size() == yp_ext_.size());

        std::copy(y_forcing.begin(), y_forcing.end(), y_ext_.begin());
        std::copy(yp_forcing.begin(), yp_forcing.end(), yp_ext_.begin());
      }

      const auto* y_system  = y_.getData();
      const auto* yp_system = yp_.getData();

      for (component_type* component : components_)
      {
        auto*       y         = component->y().getData();
        auto*       yp        = component->yp().getData();
        const auto& externals = component->getExternIndices();

        for (size_t j : externals)
        {
          IdxT local_idx = component->getNodeConnection(j);

          if (local_idx == neg1_)
          {
            y[j]  = 0.0;
            yp[j] = 0.0;
          }
          else if (local_idx < this->getInternalSize())
          {
            y[j]  = y_system[local_idx];
            yp[j] = yp_system[local_idx];
          }
          else
          {
            y[j]  = y_ext_[local_idx - this->getInternalSize()];
            yp[j] = yp_ext_[local_idx - this->getInternalSize()];
          }
        }
        component->y().setDataUpdated();
        component->yp().setDataUpdated();
      }
      return 0;
    }

    int mapGlobalToLocal()
    {
      for (auto* component : components_)
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

      for (auto* node : nodes_)
      {

        for (IdxT i = 0; i < node->size(); i++)
        {
          IdxT index = node->getNodeConnection(i);

          if (index == neg1_)
          {
            continue;
          }
          else if (internal_map_.count(index) > 0)
          {
            node->setExternalConnectionNodes(i, internal_map_.at(index));
          }
          else
          {
            node->setExternalConnectionNodes(i, external_map_.at(index));
          }
        }
      }

      return 0;
    }

    void createGlobalToInternalMap()
    {

      size_t component_internal_idx = 0;
      // First pass: Map global component indices to local subsystem indices.
      for (component_type* comp : components_)
      {
        for (IdxT i = 0; i < comp->size(); i++)
        {
          IdxT index = comp->getNodeConnection(i);

          auto extern_indices = comp->getExternIndices();

          if (index != neg1_ && !extern_indices.contains(i))
          {
            internal_map_[index] = component_internal_idx++;
          }
        }
      }

      for (node_type* node : nodes_)
      {

        for (IdxT i = 0; i < node->size(); i++)
        {
          IdxT index = node->getNodeConnection(i);

          if (index != neg1_)
          {
            internal_map_[index] = component_internal_idx++;
          }
        }
      }

      for (component_type* comp : components_)
      {
        auto extern_indices = comp->getExternIndices();
        for (IdxT j = 0; j < comp->size(); j++)
        {
          if (!extern_indices.contains(j))
          {
            continue;
          }
          IdxT index = comp->getNodeConnection(j);

          if (internal_map_.count(index) < 1 && external_map_.count(index) < 1 && index != neg1_)
          {
            external_map_[index] = component_internal_idx++;
          }
        }
      }
    }

    std::vector<IdxT>& getExternalIndices()
    {
      return external_indices_;
    }

    std::vector<ScalarT>& getExternalDataY()
    {
      return y_ext_;
    }

    std::vector<ScalarT>& getExternalDataYP()
    {
      return yp_ext_;
    }

    std::vector<ScalarT>& getExternalDataF()
    {
      return f_ext_;
    }

    void setTimeFunction(TimeFunction function)
    {
      forcing_function_ = std::move(function);
    }

    const std::unordered_map<IdxT, IdxT>& getInternalMap() const
    {
      return internal_map_;
    }

    const std::unordered_map<IdxT, IdxT>& getExternalMap() const
    {
      return external_map_;
    }

  private:
    std::unordered_map<IdxT, IdxT> internal_map_;
    std::unordered_map<IdxT, IdxT> external_map_;
    std::vector<IdxT>              external_indices_;

    std::vector<ScalarT> y_ext_;
    std::vector<ScalarT> yp_ext_;
    std::vector<ScalarT> f_ext_;

    std::optional<TimeFunction> forcing_function_;

  }; // class SubsystemModel

} // namespace GridKit