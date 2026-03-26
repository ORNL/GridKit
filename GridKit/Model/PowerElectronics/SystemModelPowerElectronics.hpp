

#pragma once

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>

#include <GridKit/Constants.hpp>
#include <GridKit/Model/PowerElectronics/Bus.hpp>
#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>
#include <GridKit/Model/PowerElectronics/CircuitNode.hpp>
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
    using CsrMatrixT     = typename CircuitComponent<ScalarT, IdxT>::CsrMatrixT;
    using component_type = CircuitComponent<ScalarT, IdxT>;
    using bus_type       = PowerElectronics::Bus<ScalarT, IdxT>;
    using node_type      = CircuitNode<ScalarT, IdxT>;

    using CircuitComponent<ScalarT, IdxT>::size_;
    using CircuitComponent<ScalarT, IdxT>::nnz_;
    using CircuitComponent<ScalarT, IdxT>::time_;
    using CircuitComponent<ScalarT, IdxT>::alpha_;
    using CircuitComponent<ScalarT, IdxT>::y_;
    using CircuitComponent<ScalarT, IdxT>::yp_;
    using CircuitComponent<ScalarT, IdxT>::f_;
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
      rel_tol_         = 1e-4;
      abs_tol_         = 1e-4;
      this->max_steps_ = 2000;
      // By default don't use the jacobian
      use_jac_         = false;
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
      rel_tol_         = rel_tol;
      abs_tol_         = abs_tol;
      this->max_steps_ = max_steps;
      // Can choose if to use jacobian
      use_jac_         = use_jac;
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
      {
        delete comp;
      }
      delete csr_jac_;
      delete[] map_to_csr_;
    }

    /**
     * @brief allocator default
     *
     * @todo this should throw an exception as no allocation without a graph is allowed.
     * Or needs to be removed from the base class
     *
     * @return int
     */
    int allocate() final
    {
      return 1;
    }

    /**
     * @brief Will check if each component has jacobian avalible. If one doesn't have it, return false.
     *
     *
     * @return true if all components have jacobian
     * @return false otherwise
     */
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

    /**
     * @brief Allocate system vectors and construct the system CSR Jacobian
     *
     * @param[in] s size of the vector (total number of unknowns)
     *
     * @post System model vectors allocated with size s
     * @post CSR Jacobian sparsity pattern is computed
     * @post COO->CSR mapping is computed
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

      // Evaluate component Jacobians to get sparsity
      distributeVectors();
      for (component_type* component : components_)
      {
        component->evaluateJacobian();
      }

      // Count the number of non-zeros
      IdxT nnz_dup = 0;
      for (const component_type* component : components_)
      {
        const IdxT* r   = component->jacobianCooRows();
        const IdxT* c   = component->jacobianCooCols();
        IdxT        nnz = component->nnz();

        for (IdxT i = 0; i < nnz; ++i)
        {
          if (component->getNodeConnection(r[i]) != neg1_ && component->getNodeConnection(c[i]) != neg1_)
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
      for (const auto& component : components_)
      {
        const IdxT*  r   = component->jacobianCooRows();
        const IdxT*  c   = component->jacobianCooCols();
        const RealT* v   = component->jacobianCooValues();
        IdxT         nnz = component->nnz();

        for (IdxT i = 0; i < nnz; ++i)
        {
          if (component->getNodeConnection(r[i]) != neg1_ && component->getNodeConnection(c[i]) != neg1_)
          {
            rows_dup[counter] = component->getNodeConnection(r[i]);
            cols_dup[counter] = component->getNodeConnection(c[i]);
            vals_dup[counter] = v[i];
            counter++;
          }
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
          if (component->getNodeConnection(j) != neg1_)
          {
            y[j]  = y_[component->getNodeConnection(j)];
            yp[j] = yp_[component->getNodeConnection(j)];
          }
          else
          {
            y[j]  = 0.0;
            yp[j] = 0.0;
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
          if (component->getNodeConnection(j) != neg1_)
          {
            f_[component->getNodeConnection(j)] += residual[j];
          }
        }
      }

      // Uncomment to print the residual in matrix market format
      // writeVectorToMatrixMarket(f_, "ScaleMicrogrid_Residual_N2_number" + std::to_string(jac_call_count_) + ".mtx", "Residual N2 number " + std::to_string(jac_call_count_));

      return 0;
    }

    /**
     * @brief Creates the system Jacobian representing \alpha dF/dy' + dF/dy
     *
     * Updates the CSR Jacobian values using the per-component mappings
     * computed during allocate().
     *
     * @return int 0 if successful, positive if there's a recoverable error, negative if unrecoverable
     */
    int evaluateJacobian() final
    {
      distributeVectors();

      // Zero out values
      RealT* vals = csr_jac_->getValues();
      for (IdxT i = 0; i < csr_jac_->getNnz(); ++i)
      {
        vals[i] = 0.0;
      }

      // Update CSR values from component Jacobians
      IdxT counter = 0;
      for (const auto& component : components_)
      {
        component->evaluateJacobian();

        const IdxT*  r   = component->jacobianCooRows();
        const IdxT*  c   = component->jacobianCooCols();
        const RealT* v   = component->jacobianCooValues();
        IdxT         nnz = component->nnz();

        for (IdxT i = 0; i < nnz; ++i)
        {
          if (component->getNodeConnection(r[i]) != neg1_ && component->getNodeConnection(c[i]) != neg1_)
          {
            vals[map_to_csr_[counter]] += v[i];
            ++counter;
          }
        }
      }

      jac_call_count_++;
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

    CsrMatrixT* getCsrJacobian() const override
    {
      return csr_jac_;
    }

    void addComponent(component_type* component)
    {
      components_.push_back(component);
    }

    void addBus(bus_type* bus)
    {
      buses_.push_back(bus);
    }

  private:
    static constexpr IdxT neg1_ = INVALID_INDEX<IdxT>;

    std::vector<component_type*> components_;
    std::vector<bus_type*>       buses_;

    IdxT*       map_to_csr_{nullptr};
    CsrMatrixT* csr_jac_{nullptr};

    int  jac_call_count_{0};
    bool use_jac_;

  }; // class PowerElectronicsModel

} // namespace GridKit
