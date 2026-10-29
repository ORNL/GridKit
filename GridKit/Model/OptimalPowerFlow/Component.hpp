/**
 * @file Component.hpp
 * @brief Base class of optimal power flow components.
 */

#pragma once

#include <cstddef>
#include <limits>
#include <numeric>
#include <string>
#include <type_traits>
#include <vector>

#include <GridKit/Constants.hpp>
#include <GridKit/LinearAlgebra/Vector/Vector.hpp>
#include <GridKit/Model/ConfigurationChecks.hpp>
#include <GridKit/Model/OptimizationEvaluator.hpp>
#include <GridKit/ScalarTraits.hpp>
#include <GridKit/Utilities/Enum.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    using Log = ::GridKit::Utilities::Logger;

    /// Limit of a variable or inequality without one
    template <typename RealT>
    inline constexpr RealT UNBOUNDED = std::numeric_limits<RealT>::infinity();

    template <typename scalar_type, typename index_type>
    class Component;

    template <typename scalar_type, typename index_type>
    class SystemModel;

    /**
     * @brief Layout of a component model: enumerations whose port counts add
     * up, and `final`, so Enzyme calls the kernels directly
     */
    template <typename ModelT, typename ScalarT, typename IdxT>
    concept ComponentKernels =
        std::derived_from<ModelT, Component<ScalarT, IdxT>>
        && std::is_final_v<ModelT>
        && Utilities::SizedEnum<typename ModelT::InternalVariablesT>
        && Utilities::SizedEnum<typename ModelT::ExternalVariablesT>
        && Utilities::SizedEnum<typename ModelT::InternalConstraintsT>
        && Utilities::SizedEnum<typename ModelT::ExternalConstraintsT>
        && ModelT::NUM_PORT_VARIABLES.size() == ModelT::NUM_PORT_CONSTRAINTS.size()
        && std::accumulate(ModelT::NUM_PORT_VARIABLES.begin(), ModelT::NUM_PORT_VARIABLES.end(), std::size_t{0})
               == Utilities::enum_size<typename ModelT::ExternalVariablesT>()
        && std::accumulate(ModelT::NUM_PORT_CONSTRAINTS.begin(), ModelT::NUM_PORT_CONSTRAINTS.end(), std::size_t{0})
               == Utilities::enum_size<typename ModelT::ExternalConstraintsT>();

    /**
     * @brief Optimal power flow component, which views its blocks of the
     * system vectors and reads and adds to those of the components at its ports
     */
    template <typename scalar_type, typename index_type>
    class Component : public Model::OptimizationEvaluator<scalar_type, index_type>
    {
      friend class SystemModel<scalar_type, index_type>;

      using EvaluatorT = Model::OptimizationEvaluator<scalar_type, index_type>;

    protected:
      using ScalarT      = scalar_type;
      using IdxT         = index_type;
      using RealT        = typename ScalarTraits<ScalarT>::RealT;
      using VectorT      = LinearAlgebra::Vector<ScalarT, IdxT>;
      using RealVectorT  = LinearAlgebra::Vector<RealT, IdxT>;
      using IndexVectorT = LinearAlgebra::Vector<IdxT, IdxT>;

    public:
      /// Numbers of variables and constraints, or the first index of each
      struct Block
      {
        IdxT variables{0};
        IdxT constraints{0};

        constexpr Block& operator+=(const Block& other)
        {
          variables   += other.variables;
          constraints += other.constraints;
          return *this;
        }
      };

      /// Blocks of `component` that a port reads and adds to
      struct Port
      {
        Block      size;
        Component* component{nullptr};
        Block      offset{};
      };

      /// Derivative entries in global indices, with room for the dense local block
      struct CooBuffer
      {
        int resize(IdxT capacity)
        {
          if (rows.resize(capacity) != 0
              || cols.resize(capacity) != 0
              || values.resize(capacity) != 0)
          {
            return 1;
          }
          return 0;
        }

        IndexVectorT rows;
        IndexVectorT cols;
        RealVectorT  values;
        IdxT         nnz{0};
      };

      virtual ~Component() = default;

      /// Every problem with the parameters
      virtual Model::ConfigurationChecks verify() const = 0;

      /// Read and add to the blocks of `component` from `offset` at `port`
      int connect(IdxT port, Component& component, const Block& offset = {});

      /// Nothing to allocate, since the system binds a component to its blocks
      int allocate() override
      {
        return 0;
      }

      /// Tag each internal constraint as an inequality or an equality
      int tagInequality() override
      {
        return 0;
      }

      /// Set the bounds and any nonzero initial value
      int initialize() override
      {
        return 0;
      }

      /// Evaluate the objective contribution into `f()`
      int evaluateObjective() override;

      /// Add the objective gradient
      int evaluateGradient() override
      {
        return 0;
      }

      /// Add the constraint contributions
      int evaluateConstraints() override;

      /// Evaluate the constraint Jacobian into `jacobianEntries()`
      int evaluateJacobian() override
      {
        return 0;
      }

      /// Evaluate the lower Lagrangian Hessian into `hessianEntries()`
      int evaluateHessian([[maybe_unused]] RealT sigma) override
      {
        return 0;
      }

      /// Objective contribution
      virtual ScalarT objective([[maybe_unused]] const ScalarT* x,
                                [[maybe_unused]] const ScalarT* external) const
      {
        return ScalarT{ZERO<RealT>};
      }

      /// Contributions `g` to the internal constraints
      virtual void internalConstraints([[maybe_unused]] const ScalarT* x,
                                       [[maybe_unused]] const ScalarT* external,
                                       [[maybe_unused]] ScalarT*       g) const
      {
      }

      /// Contributions `g` to the constraints at each port, in port order
      virtual void externalConstraints([[maybe_unused]] const ScalarT* x,
                                       [[maybe_unused]] const ScalarT* external,
                                       [[maybe_unused]] ScalarT*       g) const
      {
      }

      const std::string& id() const
      {
        return id_;
      }

      const std::vector<Port>& ports() const
      {
        return ports_;
      }

      const CooBuffer& jacobianEntries() const
      {
        return jacobian_entries_;
      }

      const CooBuffer& hessianEntries() const
      {
        return hessian_entries_;
      }

    protected:
      /// A system, which has no kernels
      Component() = default;

      /// Sizes from the enumerations of `ModelT`
      template <ComponentKernels<scalar_type, index_type> ModelT>
      Component(std::type_identity<ModelT>, std::string id);

      // Enzyme derivatives of the derived model

      template <ComponentKernels<scalar_type, index_type> ModelT>
      int evaluateGradient(const ModelT& model);

      template <ComponentKernels<scalar_type, index_type> ModelT>
      int evaluateJacobian(const ModelT& model);

      template <ComponentKernels<scalar_type, index_type> ModelT>
      int evaluateHessian(const ModelT& model, RealT sigma);

      // Blocks of the system vectors, which a component views

      using EvaluatorT::f_;
      using EvaluatorT::g_;
      using EvaluatorT::g_lower_;
      using EvaluatorT::g_upper_;
      using EvaluatorT::gradient_;
      using EvaluatorT::lambda_;
      using EvaluatorT::tag_;
      using EvaluatorT::x_;
      using EvaluatorT::x_lower_;
      using EvaluatorT::x_upper_;

    private:
      template <typename ModelT>
      static constexpr Block internalSize()
      {
        return {static_cast<IdxT>(Utilities::enum_size<typename ModelT::InternalVariablesT>()),
                static_cast<IdxT>(Utilities::enum_size<typename ModelT::InternalConstraintsT>())};
      }

      /// Internal plus external sizes of `ModelT`
      template <typename ModelT>
      static constexpr Block localSize()
      {
        return {static_cast<IdxT>(Utilities::enum_size<typename ModelT::InternalVariablesT>()
                                  + Utilities::enum_size<typename ModelT::ExternalVariablesT>()),
                static_cast<IdxT>(Utilities::enum_size<typename ModelT::InternalConstraintsT>()
                                  + Utilities::enum_size<typename ModelT::ExternalConstraintsT>())};
      }

      void gather(ScalarT* x) const;
      void gatherMultipliers(RealT* lambda) const;
      void setIndices();
      void scatterGradient(const RealT* gradient);
      void scatterConstraints(const ScalarT* g);

      static void keepLowerTriangle(CooBuffer& entries);

      std::string       id_;
      std::vector<Port> ports_;
      Block             size_;   ///< Internal sizes
      Block             offset_; ///< Global index of each internal block
      bool              bound_{false};

      VectorT      local_x_;            ///< Internal, then external variables
      VectorT      local_g_;            ///< Internal, then external contributions
      IndexVectorT variable_indices_;   ///< Global column of each local variable
      IndexVectorT constraint_indices_; ///< Global row of each local constraint contribution
      CooBuffer    jacobian_entries_;
      CooBuffer    hessian_entries_;
    };
  } // namespace OptimalPowerFlow
} // namespace GridKit
