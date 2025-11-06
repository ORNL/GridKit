

#pragma once

#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>

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

    using CircuitComponent<ScalarT, IdxT>::size_;
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
    PowerElectronicsModel()
    {
      // Set system model parameters as default
      rel_tol_                 = 1e-4;
      abs_tol_                 = 1e-4;
      this->max_steps_         = 2000;
      // By default don't use the jacobian
      use_jac_                 = false;
      jac_sparsity_row_indices = nullptr;
      jac_sparsity_col_indices = nullptr;
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
      rel_tol_                 = rel_tol;
      abs_tol_                 = abs_tol;
      this->max_steps_         = max_steps;
      // Can choose if to use jacobian
      use_jac_                 = use_jac;
      jac_sparsity_row_indices = nullptr;
      jac_sparsity_col_indices = nullptr;
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
      for (auto comp : this->components_)
        delete comp;

      if (jac_sparsity_row_indices != nullptr)
      {
        delete[] jac_sparsity_row_indices;
      }
      if (jac_sparsity_col_indices != nullptr)
      {
        delete[] jac_sparsity_col_indices;
      }
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

      // Helper struct for saving CSR sparsities of a particular component, since
      // right now they only report COO. Can be removed once components have CSR matrices.
      struct CsrSparsity
      {
        std::vector<IdxT> row_indices_;
        std::vector<IdxT> col_indices_;
      };

      // A reverse mapping from system variables -> component variables
      std::vector<std::vector<ComponentContribution>> reverse_map(size_);
      std::vector<CsrSparsity>                        component_sparsities;
      component_sparsities.reserve(components_.size());

      // Loop over all components, evaluate their jacobians, save their sparsity information,
      // and construct the reverse variable mapping.
      for (size_t comp_idx = 0; comp_idx < components_.size(); comp_idx++)
      {
        component_type* component = components_[comp_idx];
        component->evaluateJacobian();

        auto [row_indices, col_indices, _] = component->getJacobian().setDataToCSR();

        component_sparsities.push_back({
            .row_indices_ = std::move(row_indices),
            .col_indices_ = std::move(col_indices),
        });

        for (IdxT local_row = 0; local_row < component->size(); local_row++)
        {
          IdxT global_row = component->getNodeConnection(local_row);

          // Not all local variables map to a global variable
          if (global_row != neg1_)
          {
            reverse_map[global_row].push_back({.comp_idx_ = comp_idx, .local_row_idx_ = local_row});
          }
        }
      }

      // Allocate the final sparsity pattern info
      IdxT*             global_row_indices = new IdxT[size_ + 1];
      std::vector<IdxT> global_col_indices;
      global_col_indices.reserve(size_);

      // Construct sparsity pattern row-by-row. In the future, can be batched
      // into contiguous blocks of rows which can be calculated in parallel,
      // then joined sequentially
      for (size_t row = 0; row < size_; row++)
      {
        global_row_indices[row] = global_col_indices.size();

        // Collect columns from each component which has a row which contributes to this row
        for (ComponentContribution contrib : reverse_map[row])
        {
          component_type* component     = components_[contrib.comp_idx_];
          CsrSparsity&    comp_sparsity = component_sparsities[contrib.comp_idx_];

          IdxT row_start = comp_sparsity.row_indices_[contrib.local_row_idx_];
          IdxT row_end   = comp_sparsity.row_indices_[contrib.local_row_idx_ + 1];
          for (size_t local_elem_idx = row_start; local_elem_idx < row_end; local_elem_idx++)
          {
            IdxT local_col  = comp_sparsity.col_indices_[local_elem_idx];
            IdxT global_col = component->getNodeConnection(local_col);

            // Not all columns map to columns in the system jacobian
            if (global_col != neg1_)
            {
              global_col_indices.push_back(global_col);
            }
          }
        }

        // Sort the row by column indices. Since the mapping from local indices to global indices isn't monotonically increasing,
        // this is necessary.
        auto start = std::next(global_col_indices.begin(), global_row_indices[row]);
        std::sort(start, global_col_indices.end());

        // If there were multiple components which contributed columns, then we definitely need
        // to de-duplicate the columns.
        if (reverse_map[row].size() > 1)
        {
          global_col_indices.erase(std::unique(start, global_col_indices.end()), global_col_indices.end());
        }
      }
      // Final row index (beyond the end) keeps track of nnz, which is also the size of the column indices.
      global_row_indices[size_] = global_col_indices.size();

      // Delete old sparsity buffers, if they exist
      if (jac_sparsity_row_indices != nullptr)
      {
        delete[] jac_sparsity_row_indices;
      }
      if (jac_sparsity_col_indices != nullptr)
      {
        delete[] jac_sparsity_col_indices;
      }

      // Allocate new sparsity buffers
      jac_sparsity_row_indices = global_row_indices;
      jac_sparsity_col_indices = new IdxT[global_col_indices.size()];

      // Copy column indices
      std::copy(global_col_indices.begin(), global_col_indices.end(), jac_sparsity_col_indices);

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
      size_ = s;
      for (const auto& component : components_)
      {
        component->allocate();
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
     * @brief Returns a pointer to the buffer of Jacobian CSR row indices.
     * @todo Can be removed once `getJacobian()` returns a `CsrMatrix`
     */
    IdxT* getJacRowIndices()
    {
      return jac_sparsity_row_indices;
    }

    /**
     * @brief Returns a pointer to the buffer of Jacobian CSR column indices.
     * @todo Can be removed once `getJacobian()` returns a `CsrMatrix`
     */
    IdxT* getJacColIndices()
    {
      return jac_sparsity_col_indices;
    }

  private:
    static constexpr IdxT neg1_ = static_cast<IdxT>(-1);

    /**
     * @brief Keeps track of the CSR row indices of the system jacobian.
     *        `nullptr` is used to indicate an un-allocated buffer.
     * @todo  Unneeded once the jacobian is in CSR format. This can be stored
     *        in the jacobian itself.
     */
    IdxT* jac_sparsity_row_indices;
    /**
     * @brief Keeps track of the CSR col indices of the system jacobian.
     *        `nullptr` is used to indicate an un-allocated buffer.
     * @todo  Unneeded once the jacobian is in CSR format. This can be stored
     *        in the jacobian itself.
     */
    IdxT* jac_sparsity_col_indices;

    std::vector<component_type*> components_;

    int  jac_call_count_{0};
    bool use_jac_;

  }; // class PowerElectronicsModel

} // namespace GridKit
