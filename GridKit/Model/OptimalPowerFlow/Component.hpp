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
#include <GridKit/ScalarTraits.hpp>
#include <GridKit/Utilities/Enum.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    using Log = ::GridKit::Utilities::Logger;

    /// Limit of a variable or inequality constraint without one
    template <typename RealT>
    inline constexpr RealT UNBOUNDED = std::numeric_limits<RealT>::infinity();

    template <typename scalar_type, typename index_type>
    class Component;

    template <typename scalar_type, typename index_type>
    class SystemModel;

    /**
     * @brief Requirements on a component model
     *
     * A model enumerates its internal variables, the external variables it
     * reads at its ports, the internal constraints it owns, and its
     * contributions to the constraints at its ports. The counts at each port
     * add up to the external enumerations. It is `final`, so Enzyme calls its
     * kernels directly.
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
     * @brief Optimal power flow component
     *
     * `SystemModel` binds the internal blocks of each component to the system
     * vectors. Each port reads variables of another component and contributes
     * to its constraints, equalities or inequalities alike. A component
     * overrides the kernels it has; the others contribute nothing.
     */
    template <typename scalar_type, typename index_type>
    class Component
    {
      friend class SystemModel<scalar_type, index_type>;

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

      /// Variables a port reads and constraints it contributes to, in the internal blocks of `component`
      struct Port
      {
        Block      size;               ///< Numbers of variables read and constraints contributed to
        Component* component{nullptr}; ///< Owner of those variables and constraints, after `connect()`
        Block      offset{};           ///< First of them in the internal blocks of `component`
      };

      /**
       * @brief Derivative entries in global indices
       *
       * Enzyme's `sparse_store` fills the first `nnz` entries of buffers sized
       * for the dense local block. Assembly sums duplicate coordinates.
       */
      struct CooBuffer
      {
        /// Room for `capacity` entries
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

      /// Read the variables and contribute to the constraints of `component` at port `port`
      int connect(IdxT port, Component& component, const Block& offset = {});

      /// Tag the internal constraints that are inequalities; the others are equalities
      virtual int tagInequality()
      {
        return 0;
      }

      /// Set the bounds and the initial values in the internal blocks
      virtual int initialize()
      {
        return 0;
      }

      /// Evaluate the objective contribution into `f()`
      int evaluateObjective();

      /// Add the objective gradient to the system gradient
      virtual int evaluateGradient()
      {
        return 0;
      }

      /// Add the constraint contributions to the system constraints
      int evaluateConstraints();

      /// Evaluate the constraint Jacobian entries into `jacobianEntries()`
      virtual int evaluateJacobian()
      {
        return 0;
      }

      /// Evaluate the lower Lagrangian Hessian entries at the system multipliers into `hessianEntries()`
      virtual int evaluateHessian([[maybe_unused]] RealT sigma)
      {
        return 0;
      }

      /// Objective contribution over the internal variables `x` and the external variables `external`
      virtual ScalarT objective([[maybe_unused]] const ScalarT* x,
                                [[maybe_unused]] const ScalarT* external) const
      {
        return ScalarT{ZERO<RealT>};
      }

      /// Contributions `g` to the constraints this component owns
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

      /// Objective contribution at the last `evaluateObjective()`
      const ScalarT& f() const
      {
        return f_;
      }

      /// Inequality tag of each internal constraint, set by `tagInequality()`
      const std::vector<bool>& tag() const
      {
        return tag_;
      }

      /// Constraint Jacobian entries at the last `evaluateJacobian()`, in global indices
      const CooBuffer& jacobianEntries() const
      {
        return jacobian_entries_;
      }

      /// Lower Lagrangian Hessian entries at the last `evaluateHessian()`, in global indices
      const CooBuffer& hessianEntries() const
      {
        return hessian_entries_;
      }

    protected:
      /// Sizes from the enumerations that `ModelT` declares
      template <ComponentKernels<scalar_type, index_type> ModelT>
      Component(std::type_identity<ModelT>, std::string id);

      // Enzyme derivatives of `model`, the derived component

      template <ComponentKernels<scalar_type, index_type> ModelT>
      int evaluateGradient(const ModelT& model);

      template <ComponentKernels<scalar_type, index_type> ModelT>
      int evaluateJacobian(const ModelT& model);

      template <ComponentKernels<scalar_type, index_type> ModelT>
      int evaluateHessian(const ModelT& model, RealT sigma);

      /// Inequality tag of each internal constraint, which `tagInequality()` sets
      std::vector<bool> tag_;

      // Views of the internal blocks that `initialize()` sets
      VectorT     x_;
      RealVectorT x_lower_;
      RealVectorT x_upper_;
      RealVectorT g_lower_;
      RealVectorT g_upper_;

    private:
      /// Numbers of internal variables and constraints of `ModelT`
      template <typename ModelT>
      static constexpr Block internalSize()
      {
        return {static_cast<IdxT>(Utilities::enum_size<typename ModelT::InternalVariablesT>()),
                static_cast<IdxT>(Utilities::enum_size<typename ModelT::InternalConstraintsT>())};
      }

      /// Numbers of local variables and constraint contributions of `ModelT`
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
      void variableIndices(IdxT* columns) const;
      void constraintIndices(IdxT* rows) const;
      void scatterGradient(const RealT* gradient);
      void scatterConstraints(const ScalarT* g);

      static void keepLowerTriangle(CooBuffer& entries);

      std::string       id_;
      std::vector<Port> ports_;
      Block             size_;   ///< Numbers of internal variables and constraints
      Block             offset_; ///< First global index of each internal block
      bool              bound_{false};

      // Views of the internal blocks that the evaluations set, and of their multipliers
      ScalarT     f_{};
      RealVectorT gradient_;
      VectorT     g_;
      RealVectorT lambda_;

      VectorT   local_x_; ///< Local variables: internal, then external
      VectorT   local_g_; ///< Local constraint contributions: internal, then external
      CooBuffer jacobian_entries_;
      CooBuffer hessian_entries_;
    };
  } // namespace OptimalPowerFlow
} // namespace GridKit
