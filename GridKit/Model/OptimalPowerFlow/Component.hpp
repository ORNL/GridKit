/**
 * @file Component.hpp
 * @brief Base class of optimal power flow components.
 */

#pragma once

#include <array>
#include <concepts>
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
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    using Log = ::GridKit::Utilities::Logger;

    /// Limit of a variable or inequality constraint function without one
    template <typename RealT>
    inline constexpr RealT UNBOUNDED = std::numeric_limits<RealT>::infinity();

    /// Sum of the counts at each port
    template <typename IdxT, std::size_t N>
    constexpr IdxT total(const std::array<IdxT, N>& counts)
    {
      return std::accumulate(counts.begin(), counts.end(), IdxT{0});
    }

    template <typename scalar_type, typename index_type>
    class Component;

    template <typename scalar_type, typename index_type>
    class SystemModel;

    /**
     * @brief Requirements on a component
     *
     * A component declares its ports, its internal variables and rows, and at
     * each port the variables it reads and the equality rows it contributes
     * to. It is `final`, so Enzyme calls its kernels directly. Local variables
     * and equality contributions are internal first, then those of each port.
     */
    template <typename ModelT, typename ScalarT, typename IdxT>
    concept ComponentKernels =
        std::derived_from<ModelT, Component<ScalarT, IdxT>>
        && std::is_final_v<ModelT>
        && ModelT::NUM_PORT_VARIABLES.size() == ModelT::NUM_PORTS
        && ModelT::NUM_PORT_EQUALITIES.size() == ModelT::NUM_PORTS
        && requires {
             ModelT::NUM_INTERNAL_VARIABLES;
             ModelT::NUM_INTERNAL_EQUALITIES;
             ModelT::NUM_INTERNAL_INEQUALITIES;
           };

    /**
     * @brief Optimal power flow component
     *
     * `SystemModel` binds the internal blocks of each component to the system
     * vectors. Each port reads variables of another component and contributes
     * to its equality rows. A component overrides the kernels it has; the
     * others contribute nothing.
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
      /// Numbers of variables, equality rows, and inequality rows, or the first global index of each
      struct Block
      {
        IdxT variables{0};
        IdxT equalities{0};
        IdxT inequalities{0};

        constexpr Block& operator+=(const Block& other)
        {
          variables    += other.variables;
          equalities   += other.equalities;
          inequalities += other.inequalities;
          return *this;
        }
      };

      /// Variables a port reads and equality rows it contributes to, in the internal blocks of `component`
      struct Port
      {
        Block      size;               ///< Numbers of variables read and equality rows contributed to
        Component* component{nullptr}; ///< Owner of those variables and rows, after `connect()`
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

      /// Components set their limits in their views
      virtual int allocate()
      {
        return 0;
      }

      /// Components set their starting point in their views
      virtual int initialize()
      {
        return 0;
      }

      int evaluateObjective();
      int evaluateEqualities();
      int evaluateInequalities();

      virtual int evaluateGradient()
      {
        return 0;
      }

      virtual int evaluateEqualityJacobian()
      {
        return 0;
      }

      virtual int evaluateInequalityJacobian()
      {
        return 0;
      }

      virtual int evaluateHessian([[maybe_unused]] RealT        sigma,
                                  [[maybe_unused]] const RealT* lambda_c,
                                  [[maybe_unused]] const RealT* lambda_d)
      {
        return 0;
      }

      /// Objective contribution over the local variables `x`
      virtual ScalarT objective([[maybe_unused]] const ScalarT* x) const
      {
        return ScalarT{ZERO<RealT>};
      }

      /// Contributions `c` to the equality rows this component owns
      virtual void internalEqualities([[maybe_unused]] const ScalarT* x, [[maybe_unused]] ScalarT* c) const
      {
      }

      /// Inequality rows `d` this component owns, over the local variables `x`
      virtual void internalInequalities([[maybe_unused]] const ScalarT* x, [[maybe_unused]] ScalarT* d) const
      {
      }

      /// Contributions `c` to the equality rows at each port, in port order
      virtual void externalEqualities([[maybe_unused]] const ScalarT* x, [[maybe_unused]] ScalarT* c) const
      {
      }

      /// Read the variables and contribute to the equality rows of `component` at port `port`
      int connect(IdxT port, Component& component, const Block& offset = {});

      /// Objective contribution at the last `evaluateObjective()`
      const ScalarT& f() const
      {
        return f_;
      }

      const std::string& id() const
      {
        return id_;
      }

      const std::vector<Port>& ports() const
      {
        return ports_;
      }

      /// Equality Jacobian entries at the last `evaluateEqualityJacobian()`, in global indices
      const CooBuffer& equalityEntries() const
      {
        return equality_entries_;
      }

      /// Inequality Jacobian entries at the last `evaluateInequalityJacobian()`, in global indices
      const CooBuffer& inequalityEntries() const
      {
        return inequality_entries_;
      }

      /// Lower Lagrangian Hessian entries at the last `evaluateHessian()`, in global indices
      const CooBuffer& hessianEntries() const
      {
        return hessian_entries_;
      }

    protected:
      /// Sizes from the counts that `ModelT` declares
      template <ComponentKernels<scalar_type, index_type> ModelT>
      Component(std::type_identity<ModelT>, std::string id);

      // Enzyme derivatives of `model`, the derived component
      template <ComponentKernels<scalar_type, index_type> ModelT>
      int evaluateGradient(const ModelT& model);

      template <ComponentKernels<scalar_type, index_type> ModelT>
      int evaluateEqualityJacobian(const ModelT& model);

      template <ComponentKernels<scalar_type, index_type> ModelT>
      int evaluateInequalityJacobian(const ModelT& model);

      template <ComponentKernels<scalar_type, index_type> ModelT>
      int evaluateHessian(const ModelT& model, RealT sigma, const RealT* lambda_c, const RealT* lambda_d);

      ScalarT     f_{};
      VectorT     x_;
      RealVectorT x_lower_;
      RealVectorT x_upper_;
      RealVectorT gradient_;
      VectorT     c_;
      VectorT     d_;
      RealVectorT d_lower_;
      RealVectorT d_upper_;

    private:
      /// Numbers of local variables, equality contributions, and inequality rows of `ModelT`
      template <typename ModelT>
      static constexpr Block localSize()
      {
        return {ModelT::NUM_INTERNAL_VARIABLES + total(ModelT::NUM_PORT_VARIABLES),
                ModelT::NUM_INTERNAL_EQUALITIES + total(ModelT::NUM_PORT_EQUALITIES),
                ModelT::NUM_INTERNAL_INEQUALITIES};
      }

      void gather(ScalarT* x) const;
      void globalIndices(IdxT Block::* block, IdxT* indices) const;

      template <typename ValueT, typename MemberT>
      void scatter(IdxT Block::* block, MemberT Component::* vector, const ValueT* local);

      template <typename DriverT, typename ModelT>
      int evaluateJacobian(const ModelT& model, const IdxT* rows, CooBuffer& entries);

      template <typename Keep>
      static void filter(CooBuffer& entries, Keep keep);

      std::string       id_;
      std::vector<Port> ports_;
      Block             size_;   ///< Numbers of internal variables and rows
      Block             offset_; ///< First global index of each internal block
      bool              bound_{false};
      VectorT           local_x_; ///< Local variables
      VectorT           local_c_; ///< Local equality contributions
      CooBuffer         equality_entries_;
      CooBuffer         inequality_entries_;
      CooBuffer         hessian_entries_;
    };
  } // namespace OptimalPowerFlow
} // namespace GridKit
