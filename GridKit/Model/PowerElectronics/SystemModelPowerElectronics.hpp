

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
        delete comp;
    }

    /**
     * @brief allocator default
     *
     * @todo this should throw an exception as no allocation without a graph is allowed.
     * Or needs to be removed from the base class
     *
     * @return int
     */
    int allocate()
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

  private:
    static constexpr IdxT neg1_ = static_cast<IdxT>(-1);

    std::vector<component_type*> components_;

    int  jac_call_count_{0};
    bool use_jac_;

  }; // class PowerElectronicsModel

} // namespace GridKit
