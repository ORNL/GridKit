

#pragma once

#include <algorithm>
#include <cassert>
#include <forward_list>
#include <iomanip>
#include <iostream>
#include <vector>

#include <GridKit/LinearAlgebra/SparseMatrix/CsrMatrix.hpp>
#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>
#include <GridKit/Model/PowerElectronics/CircuitGraph.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{
  /**
   * Writes a vector to a file in Matrix Market format
   *
   * @param vec The vector to write
   * @param filename The name of the output file
   * @param header Additional header information/comments
   * @return true if the write was successful, false otherwise
   */
  template <typename T>
  void writeVectorToMatrixMarket(const std::vector<T>& vec, const std::string& filename, const std::string& header)
  {
    std::ofstream outFile(filename);

    if (!outFile.is_open())
    {
      std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
      return;
    }

    // Uncomment to write Matrix Market header
    // outFile << "%%MatrixMarket vector array real general" << std::endl;

    // Write additional header information as comments
    if (!header.empty())
    {
      outFile << "% " << header << std::endl;
    }

    // Write the vector size
    outFile << vec.size() << std::endl;

    // Write the vector elements
    outFile << std::scientific << std::setprecision(16);
    for (const auto& val : vec)
    {
      outFile << val << std::endl;
    }

    outFile.close();
    return;
  }

  template <class ScalarT, typename IdxT>
  class PowerElectronicsModel : public CircuitComponent<ScalarT, IdxT>
  {
    using RealT          = typename CircuitComponent<ScalarT, IdxT>::RealT;
    using component_type = CircuitComponent<ScalarT, IdxT>;
    using MatrixT        = CircuitComponent<ScalarT, IdxT>::MatrixT;

    using CircuitComponent<ScalarT, IdxT>::size_;
    using CircuitComponent<ScalarT, IdxT>::n_intern_;
    using CircuitComponent<ScalarT, IdxT>::n_extern_;
    using CircuitComponent<ScalarT, IdxT>::nnz_;
    using CircuitComponent<ScalarT, IdxT>::time_;
    using CircuitComponent<ScalarT, IdxT>::alpha_;
    using CircuitComponent<ScalarT, IdxT>::y_;
    using CircuitComponent<ScalarT, IdxT>::yp_;
    using CircuitComponent<ScalarT, IdxT>::f_;
    using CircuitComponent<ScalarT, IdxT>::jac_;
    using CircuitComponent<ScalarT, IdxT>::rel_tol_;
    using CircuitComponent<ScalarT, IdxT>::abs_tol_;
    using CircuitComponent<ScalarT, IdxT>::max_steps_;

  public:
    /**
     * @brief Default constructor for the system model
     *
     * @post System model parameters set as default
     */
    PowerElectronicsModel()
    {
      // Set system model parameters as default
      rel_tol_   = 1e-4;
      abs_tol_   = 1e-4;
      max_steps_ = 2000;
      // By default don't use the jacobian
      use_jac_   = false;
    }

    /**
     * @brief Constructor for the system model
     *
     * @param[in] rel_tol Relative tolerance for the system model
     * @param[in] abs_tol Absolute tolerance for the system model
     * @param[in] use_jac Boolean to choose if to use jacobian
     * @param[in] max_steps Maximum number of steps for the system model
     *
     * @post System model parameters set as input
     */
    PowerElectronicsModel(double rel_tol   = 1e-4,
                          double abs_tol   = 1e-4,
                          bool   use_jac   = false,
                          IdxT   max_steps = 2000)
    {
      // Set system model tolerances from input
      rel_tol_   = rel_tol;
      abs_tol_   = abs_tol;
      max_steps_ = max_steps;
      // Can choose if to use jacobian
      use_jac_   = use_jac;
    }

    /**
     * @brief Destructor for the system model
     *
     * @pre System components are allocated
     *
     * @post System components are deallocated
     *
     */
    virtual ~PowerElectronicsModel()
    {
      for (auto comp : components_)
        delete comp;
    }

    /**
     * @brief Allocates and constructs sparsity patterns for system jacobian.
     *        To do this, sparsity patterns of components are needed, so
     *        each component's jacobian is evaluated, and the sparsity pattern
     *        of that component's jacobian is expected to match the sparsity pattern
     *        for the rest of the simulation.
     * @pre   `size_` must be correctly initialized
     */
    int allocate()
    {
      // Helper struct for identifying a particular component's local system variable
      struct ComponentContribution
      {
        // The index of the component in `components_`.
        size_t comp_idx_;
        // The local system variable index
        size_t local_row_idx_;
      };

      // A reverse mapping from external system variables -> component variables
      std::forward_list<ComponentContribution>* reverse_extern_map = new std::forward_list<ComponentContribution>[n_extern_];
      size_t                                    component_nnz      = 0;

      // Loop over all components, evaluate their jacobians, save their sparsity information,
      // and construct the reverse variable mapping.
      for (size_t comp_idx = 0; comp_idx < components_.size(); comp_idx++)
      {
        component_type* component = components_[comp_idx];
        component->evaluateJacobian();

        MatrixT& comp_jacobian = component->getJacobian();

        for (IdxT local_external_row : component->getExternIndices())
        {
          IdxT global_row = component->getNodeConnection(local_external_row);

          // Not all local variables map to a global variable
          if (global_row != neg1_)
          {
            // global_row >= n_intern_ (see allocate(s))
            reverse_extern_map[global_row - n_intern_].push_front({.comp_idx_ = comp_idx, .local_row_idx_ = local_external_row});
          }
        }

        component_nnz += comp_jacobian.nnz();
      }

      // Allocate the final sparsity pattern info
      IdxT* global_row_indices = new IdxT[size_ + 1];
      IdxT* global_col_indices = new IdxT[component_nnz]; // Use component_nnz as an upper bound on nnz
      global_row_indices[0]    = 0;

      // Construct Jacobian sparsity pattern

      // Start with internal rows, which must be before external rows, are grouped by component,
      // strictly increasing wrt component internals, and must be sorted by component (see allocate(s))
      size_t curr_internal_row = 0;
      for (size_t comp_idx = 0; comp_idx < components_.size(); comp_idx++)
      {
        component_type*                                                         component      = components_[comp_idx];
        std::set<size_t>                                                        comp_externals = component->getExternIndices();
        std::tuple<std::vector<IdxT>&, std::vector<IdxT>&, std::vector<RealT>&> entries        = component->getJacobian().getEntries();
        const std::vector<IdxT>&                                                rows           = std::get<0>(entries);
        const std::vector<IdxT>&                                                columns        = std::get<1>(entries);

        if (rows.empty())
          continue;

        // Loop through all elements of the component jacobian, adding them to the system jacobian if needed
        for (size_t elem_idx = 0; elem_idx < rows.size(); elem_idx++)
        {
          IdxT local_row  = rows[elem_idx];
          IdxT global_row = component->getNodeConnection(local_row);

          // Skip variables which aren't system variables
          if (comp_externals.contains(local_row) || global_row == neg1_)
            continue;

          if (global_row > curr_internal_row)
          {
            curr_internal_row++;
            global_row_indices[curr_internal_row + 1] = global_row_indices[curr_internal_row];
          }

          assert(global_row == curr_internal_row);

          IdxT local_col_idx  = columns[elem_idx];
          IdxT global_col_idx = component->getNodeConnection(local_col_idx);

          if (global_col_idx == neg1_)
            continue;

          global_col_indices[global_row_indices[curr_internal_row + 1]] = global_col_idx;
          global_row_indices[curr_internal_row + 1]++;
        }
      }

      // One last row after
      curr_internal_row++;
      global_row_indices[curr_internal_row + 1] = global_row_indices[curr_internal_row];

      // Need to sort internal rows because even though the mapping from component internal to system internal
      // is monotonic, the mapping from component external to system external may not be, and internal rows
      // may contain external columns.
      for (size_t row_idx = 0; row_idx < curr_internal_row; row_idx++)
      {
        IdxT* global_row_start = global_col_indices + global_row_indices[row_idx];
        IdxT* global_row_end   = global_col_indices + global_row_indices[row_idx + 1];
        std::sort(global_row_start, global_row_end);
      }

      assert(curr_internal_row == n_intern_);

      // Then move on to external rows, which come after internal rows
      for (size_t row = n_intern_; row < size_; row++)
      {
        global_row_indices[row + 1] = global_row_indices[row];

        // Collect columns from each component which has a row which contributes to this row
        for (ComponentContribution contrib : reverse_extern_map[row - n_intern_])
        {
          component_type*                                                         component = components_[contrib.comp_idx_];
          std::tuple<std::vector<IdxT>&, std::vector<IdxT>&, std::vector<RealT>&> entries   = component->getJacobian().getEntries();
          const std::vector<IdxT>&                                                rows      = std::get<0>(entries);
          const std::vector<IdxT>&                                                columns   = std::get<1>(entries);

          typename std::vector<IdxT>::const_iterator row_start = std::ranges::lower_bound(rows, contrib.local_row_idx_);

          // It can happen where external contributions are only constants, and do not appear in the jacobian.
          // If that is the case, we won't be able to find local_row_idx_ and must skip this contribution
          if (row_start == rows.end() || *row_start != contrib.local_row_idx_)
            continue;

          typename std::vector<IdxT>::const_iterator row_end = std::upper_bound(row_start, rows.end(), contrib.local_row_idx_);

          for (size_t local_elem_idx = std::distance(rows.begin(), row_start);
               local_elem_idx < static_cast<size_t>(std::distance(rows.begin(), row_end));
               local_elem_idx++)
          {
            IdxT local_col  = columns[local_elem_idx];
            IdxT global_col = component->getNodeConnection(local_col);

            // Not all columns map to columns in the system jacobian
            if (global_col != neg1_)
            {
              global_col_indices[global_row_indices[row + 1]] = global_col;
              global_row_indices[row + 1]++;
            }
          }
        }

        // Sort the row by column indices. Since the mapping from local indices to global indices isn't monotonically increasing,
        // this is necessary.
        IdxT* start = global_col_indices + global_row_indices[row];
        IdxT* end   = global_col_indices + global_row_indices[row + 1];
        std::sort(start, end);

        // De-duplicate the columns
        IdxT* new_end               = std::unique(start, end);
        global_row_indices[row + 1] = global_row_indices[row] + static_cast<IdxT>(std::distance(start, new_end));
      }
      // Allocate new sparsity buffers
      IdxT nnz = global_row_indices[size_];

      csr_jac_.resize(size_, size_);
      csr_jac_.setNnz(nnz);
      csr_jac_.allocateMatrixData(LinearAlgebra::memory::HOST);

      // Copy column indices
      std::copy(global_col_indices, global_col_indices + nnz, csr_jac_.getColData());
      std::copy(global_row_indices, global_row_indices + size_ + 1, csr_jac_.getRowData());

      delete[] reverse_extern_map;
      delete[] global_row_indices;
      delete[] global_col_indices;

      return 1;
    }

    /**
     * @brief Will check if each component has jacobian avalible. If one doesn't have it, return false.
     *
     *
     * @return true if all components have jacobian
     * @return false otherwise
     */
    bool hasJacobian()
    {
      if (!use_jac_)
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

    /**
     * @brief Allocate the vector data with size amount
     * @todo Add capability to go through component model connection to get the size of the actual vector
     *
     * @param[in] s size of the vector
     *
     * @post System model vectors allocated with size s
     *
     * @return int 0 if successful, positive if there's a recoverable error, negative if unrecoverable
     */
    int allocate(IdxT s)
    {
      // Allocate all components
      size_     = s;
      n_intern_ = 0;
      for (const auto& component : components_)
      {
        component->allocate();

        // Count up the amount of internal variables which get mapped to system variables
        std::set<IdxT> extern_indices = component->getExternIndices();
        for (IdxT comp_var_idx = 0; comp_var_idx < component->size(); comp_var_idx++)
        {
          IdxT sys_var_idx = component->getNodeConnection(comp_var_idx);
          if (!extern_indices.contains(comp_var_idx) && sys_var_idx != neg1_)
          {
            n_intern_++;
          }
        }
      }
      n_extern_ = size_ - n_intern_;

      // Ensure that all component locals are mapped to system locals
      // and all component externals are mapped to system externals.
      // System locals are stored first, in 0..n_intern_ and externals
      // are stored after, in n_intern_..
      // Additionally, ensure that components locals are grouped in the system vectors - no other
      // variables are between locals from a single component. As well, ensure that components are
      // sorted by these groupings, so the first component is the first block and so on.
      for (size_t comp_idx = 0; comp_idx < components_.size(); comp_idx++)
      {
        component_type* component      = components_[comp_idx];
        std::set<IdxT>  extern_indices = component->getExternIndices();

        // Whether or not we've seen a local variable yet
        bool has_seen_local = false;
        // The system index of the last local we've seen
        IdxT last_local_sys;

        for (IdxT comp_var_idx = 0; comp_var_idx < component->size(); comp_var_idx++)
        {
          IdxT sys_var_idx = component->getNodeConnection(comp_var_idx);

          // Ignore local variables which aren't mapped to sytem variables
          if (sys_var_idx == neg1_)
            continue;

          if (extern_indices.contains(comp_var_idx))
          {
            assert(sys_var_idx >= n_intern_);
          }
          else
          {
            assert(sys_var_idx < n_intern_);

            // Make sure that all of the locals for a particular component are grouped
            // together in the sytem vector, and have increasing indices.
            if (has_seen_local)
            {
              assert(sys_var_idx == last_local_sys + 1);
            }
            else if (comp_idx > 0)
            {
              // Ensure increasing indices in components - so no need to sort components
              assert(sys_var_idx > last_local_sys);
            }

            has_seen_local = true;
            last_local_sys = sys_var_idx;
          }
        }
      }

      // Allocate global vectors
      y_.resize(size_);
      yp_.resize(size_);
      f_.resize(size_);

      allocate();

      return 0;
    }

    /**
     * @brief Set intial y and y' of each component
     *
     * @return int 0 if successful, positive if there's a recoverable error, negative if unrecoverable
     */
    int initialize()
    {
      // Initialize components
      for (const auto& component : components_)
      {
        component->initialize();
      }
      distributeVectors();

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
        for (IdxT j = 0; j < component->size(); ++j)
        {
          if (component->getNodeConnection(j) != neg1_)
          {
            component->y()[j]  = y_[component->getNodeConnection(j)];
            component->yp()[j] = yp_[component->getNodeConnection(j)];
          }
          else
          {
            component->y()[j]  = 0.0;
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
     * @return int 0 if successful, positive if there's a recoverable error, negative if unrecoverable
     */
    int evaluateResidual()
    {
      for (IdxT i = 0; i < f_.size(); i++)
      {
        f_[i] = 0.0;
      }

      distributeVectors();

      // Update system residual vector

      for (const auto& component : components_)
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

      // Uncomment to print the residual in matrix market format
      // writeVectorToMatrixMarket(f_, "ScaleMicrogrid_Residual_N2_number" + std::to_string(jac_call_count_) + ".mtx", "Residual N2 number " + std::to_string(jac_call_count_));

      return 0;
    }

    /**
     * @brief Creates the Sparse COO Jacobian representing  \alpha dF/dy' + dF/dy
     *
     * @return int 0 if successful, positive if there's a recoverable error, negative if unrecoverable
     */
    int evaluateJacobian()
    {
      jac_.zeroMatrix();
      distributeVectors();

      // Evaluate component jacs
      for (const auto& component : components_)
      {
        component->evaluateJacobian();

        // get references to local jacobian
        std::tuple<std::vector<IdxT>&, std::vector<IdxT>&, std::vector<RealT>&> tpm = component->getJacobian().getEntries();
        const auto& [r, c, v]                                                       = tpm;

        // Create copies of data to handle groundings
        std::vector<IdxT>  rgr;
        std::vector<IdxT>  cgr;
        std::vector<RealT> vgr;
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
        // elementwise jac_(rgr, cgr) += vgr
        jac_.axpy(1.0, rgr, cgr, vgr);
      }

      // jac_.printMatrixMarket("ScaleMicrogrid_Jacobian_N2_number" + std::to_string(jac_call_count_) + ".mtx", "Jacobian N2 number " + std::to_string(jac_call_count_));
      jac_call_count_++;

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
    void updateTime(RealT t, RealT a)
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

    void addComponent(component_type* component)
    {
      components_.push_back(component);
    }

    /**
     * @brief Returns a reference to the CSR Jacobian.
     * @todo Can be removed once `getJacobian()` returns a `CsrMatrix`
     */
    LinearAlgebra::CsrMatrix<RealT, IdxT>& getCsrJac()
    {
      return csr_jac_;
    }

  private:
    static constexpr IdxT neg1_ = static_cast<IdxT>(-1);

    LinearAlgebra::CsrMatrix<RealT, IdxT> csr_jac_;

    std::vector<component_type*> components_;

    int  jac_call_count_{0};
    bool use_jac_;

  }; // class PowerElectronicsModel

} // namespace GridKit
