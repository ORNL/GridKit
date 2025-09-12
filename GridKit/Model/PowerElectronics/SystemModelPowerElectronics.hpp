

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

    using CircuitComponent<ScalarT, IdxT>::size;

    using typename Model::Evaluator<ScalarT, IdxT>::CsrJacobian;

    /**
     * @brief A plan for constructing the system jacobian. Consists of individual plans
     * for constructing each row of the jacobian, which must be executed in order to construct
     * the overall matrix.
     */
    struct JacobianAssemblyPlan
    {
      /**
       * @brief A plan for constructing a row of the system jacobian which corresponds to an internal equation
       */
      struct InternalRowPlan
      {
        /**
         * @brief The index of the component that this equation belongs to
         */
        size_t            component_idx_;
        /**
         * @brief The index of the row of the component's jacobian that this row corresponds to.
         */
        size_t            row_idx_;
        /**
         * @brief The expected number of nonzero elements in the row of the component's jacobian.
         */
        IdxT              row_nnz_;
        /**
         * @brief The column indices of the elements in this row of the system jacobian
         */
        std::vector<IdxT> col_indices_;
      };

      /**
       * @brief A plan for constructing a row of the system jacobian which corresponds to an external equation
       */
      struct ExternalRowPlan
      {
        /**
         * @brief A single element in the row. An element is the sum of the elements of (potentially) several component jacobians
         */
        struct SystemElement
        {
          /**
           * @brief A single element in a component's row. Represents a summand in the calculation of a single system jacobian element
           * in an external row.
           */
          struct ComponentElement
          {
            /**
             * @brief The index of the component that this element belongs to.
             */
            size_t component_idx_;
            /**
             * @brief The index of the element inside the component jacobian's `CsrMatrix::values()` and `CsrMatrix::column_indices()` lists.
             */
            size_t element_idx_;
            /**
             * @brief The expected column of the element inside the component jacobian.
             */
            IdxT   column_idx_;
          };

          /**
           * @brief A list of all elements in component jacobians which sum to this element in the system jacobian.
           */
          std::vector<ComponentElement> component_elements_;
          /**
           * @brief The column of this element in the system jacobian.
           */
          IdxT                          column_;
        };

        /**
         * @brief A list of all elements on this row of the system jacobian.
         */
        std::vector<SystemElement> elements_;
      };

      using RowPlan = std::variant<InternalRowPlan, ExternalRowPlan>;

      /**
       * @brief A plan for assembling each row
       */
      std::vector<RowPlan> row_plans_;

      /**
       * @brief Number of components which are part of the system.
       */
      size_t num_components_;

      /**
       * @brief The dimension of the system - i.e. how many variables/residual there are
       */
      size_t system_size_;
    };

    struct ComponentCsrJacobianView
    {
      std::vector<std::variant<CsrJacobian*, CsrJacobian>> component_jacobians_;

      const CsrJacobian& operator[](size_t component_idx) const
      {
        auto& either = component_jacobians_[component_idx];
        try
        {
          return *std::get<CsrJacobian*>(either);
        }
        catch (...)
        {
          return std::get<CsrJacobian>(either);
        }
      }
    };

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

    bool hasCsrJacobian() override
    {
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
      for (auto component : components_)
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

    int evaluateCsrJacobian() override
    {
      using ExternalRowPlan = JacobianAssemblyPlan::ExternalRowPlan;
      using InternalRowPlan = JacobianAssemblyPlan::InternalRowPlan;

      auto apply_internal_row_plan = [](const InternalRowPlan& plan, const ComponentCsrJacobianView& jac_view, auto& builder)
      {
        const CsrJacobian& component_jac = jac_view[plan.component_idx_];

        IdxT row_start = component_jac.rowIndices()[plan.row_idx_];
        IdxT row_end   = component_jac.rowIndices()[plan.row_idx_ + 1];

        assert(component_jac.rowIndices().size() >= plan.row_idx_ + 2);
        assert(row_end - row_start == plan.row_nnz_);

        for (IdxT i = row_start; i < row_end; i++)
        {
          builder.elem(plan.col_indices_[i - row_start], component_jac.values()[i]);
        }
      };

      auto apply_external_row_plan = [](const ExternalRowPlan& plan, const ComponentCsrJacobianView& jac_view, auto& builder)
      {
        for (auto& element : plan.elements_)
        {
          ScalarT sum = 0;

          for (auto& element : element.component_elements_)
          {
            const CsrJacobian& component_jac = jac_view[element.component_idx_];

            assert(component_jac.colIndices()[element.element_idx_] == element.column_idx_);

            sum += component_jac.values()[element.element_idx_];
          }

          builder.elem(element.column_, sum);
        }
      };

      auto apply_jacobian_assembly_plan = [&](const JacobianAssemblyPlan& plan, const ComponentCsrJacobianView& jac_view, auto builder)
      {
        assert(size() == plan.system_size_);
        assert(components_.size() == plan.num_components_);

        for (size_t row = 0; row < plan.row_plans_.size(); row++)
        {
          builder.row(row);

          auto row_plan = plan.row_plans_[row];
          try
          {
            apply_internal_row_plan(std::get<InternalRowPlan>(row_plan), jac_view, builder);
          }
          catch (...)
          {
            apply_external_row_plan(std::get<ExternalRowPlan>(row_plan), jac_view, builder);
          }
        }

        csr_jacobian_ = std::move(builder);
      };

      auto component_jac_view = createComponentCsrJacobianView();

      using CsrBuilder = LinearAlgebra::CsrBuilder<ScalarT, IdxT, true, false>;

      if (jacobian_assembly_plan_)
      {
        apply_jacobian_assembly_plan(*jacobian_assembly_plan_, component_jac_view, CsrBuilder::fromTemplate(std::move(csr_jacobian_)));
      }

      if (!jacobian_assembly_plan_)
      {
        createJacobianAssemblyPlan(component_jac_view);
        apply_jacobian_assembly_plan(*jacobian_assembly_plan_, component_jac_view, CsrBuilder::fromEmpty(size(), size()));
      }

      return 0;
    }

    CsrJacobian& getCsrJacobian() override
    {
      return csr_jacobian_;
    }

    const CsrJacobian& getCsrJacobian() const override
    {
      return csr_jacobian_;
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

    std::vector<component_type*>        components_;
    std::optional<JacobianAssemblyPlan> jacobian_assembly_plan_;
    CsrJacobian                         csr_jacobian_{0, 0};

    int  jac_call_count_{0};
    bool use_jac_;

    ComponentCsrJacobianView createComponentCsrJacobianView()
    {
      ComponentCsrJacobianView view;

      view.component_jacobians_.reserve(components_.size());

      for (auto component : components_)
      {
        if (component->hasCsrJacobian())
        {
          component->evaluateCsrJacobian();
          view.component_jacobians_.push_back(&component->getCsrJacobian());
        }
        else
        {
          component->evaluateJacobian();
          view.component_jacobians_.push_back(CsrJacobian::fromCOO(component->getJacobian()));
        }
      }

      return view;
    }

    void createJacobianAssemblyPlan(const ComponentCsrJacobianView& component_jacs)
    {
      JacobianAssemblyPlan plan;

      plan.num_components_ = components_.size();
      plan.system_size_    = size();

      auto inverse_map = inverseComponentConnectionMap();

      for (size_t row = 0; row < plan.system_size_; row++)
      {
        auto component_contributions = inverse_map[row];
        if (component_contributions.size() > 1)
        {
          plan.row_plans_[row] = createExternalRowPlan(component_jacs, row, component_contributions);
        }
        else
        {
          auto [comp_idx, local_idx] = component_contributions.front();
          plan.row_plans_[row]       = createInternalRowPlan(component_jacs[comp_idx], row, comp_idx, local_idx);
        }
      }

      jacobian_assembly_plan_ = plan;
    }

    /**
     * @brief Returns a vector which indicates which system equations are external (true)
     */
    std::vector<std::vector<std::tuple<size_t, IdxT>>> inverseComponentConnectionMap()
    {
      std::vector<std::vector<std::tuple<size_t, IdxT>>> map(size());

      for (size_t comp_idx = 0; comp_idx < components_.size(); comp_idx++)
      {
        auto comp = components_[comp_idx];
        for (IdxT local_idx = 0; local_idx < comp->size(); local_idx++)
        {
          size_t global_idx = static_cast<size_t>(comp->getNodeConnection(local_idx));

          map[global_idx].push_back({comp_idx, local_idx});
        }
      }

      return map;
    }

    JacobianAssemblyPlan::ExternalRowPlan createExternalRowPlan(
        const ComponentCsrJacobianView&              component_jacs,
        size_t                                       row,
        const std::vector<std::tuple<size_t, IdxT>>& component_contributions)
    {
      std::vector<std::vector<typename JacobianAssemblyPlan::ExternalRowPlan::SystemElement::ComponentElement>> component_elements(size());
      size_t                                                                                                    nnz = 0;

      for (const auto& contribution : component_contributions)
      {
        auto [comp_idx, local_idx] = contribution;
        const auto& comp_jac       = component_jacs[comp_idx];

        for (size_t elem_idx = comp_jac.rowIndices()[local_idx]; elem_idx < comp_jac.rowIndices()[local_idx + 1]; elem_idx++)
        {
          IdxT   local_col_idx  = comp_jac.colIndices()[elem_idx];
          size_t global_col_idx = components_[comp_idx]->getNodeConnection(local_col_idx);

          if (component_elements[global_col_idx].empty())
          {
            nnz++;
          }

          component_elements[global_col_idx].push_back({.component_idx_ = comp_idx,
                                                        .element_idx_   = elem_idx,
                                                        .column_idx_    = local_col_idx});
        }
      }

      typename JacobianAssemblyPlan::ExternalRowPlan rowPlan;
      rowPlan.elements_.reserve(nnz);

      for (size_t col_idx = 0; col_idx < component_elements.size(); col_idx++)
      {
        if (!component_elements[col_idx].empty())
        {
          rowPlan.elements_.push_back({
              .component_elements_ = std::move(component_elements[col_idx]),
              .column_             = col_idx,
          });
        }
      }

      return rowPlan;
    }

    JacobianAssemblyPlan::InternalRowPlan createInternalRowPlan(
        const CsrJacobian& component_jac,
        size_t             row,
        size_t             comp_idx,
        IdxT               local_idx) const
    {
      IdxT        row_idx      = component_jac.rowIndices()[row];
      IdxT        next_row_idx = component_jac.rowIndices()[row + 1];
      std::vector col_indices(std::next(component_jac.colIndices().begin(), row_idx), std::next(component_jac.colIndices().begin(), next_row_idx));

      return {
          .component_idx_ = comp_idx,
          .row_idx_       = local_idx,
          .row_nnz_       = next_row_idx - row_idx,
          .col_indices_   = col_indices,
      };
    }

  }; // class PowerElectronicsModel

} // namespace GridKit
