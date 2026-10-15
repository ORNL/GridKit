/**
 * @file Component.hpp
 * @brief Base class of optimal power flow components.
 */

#pragma once

#include <limits>
#include <string>
#include <utility>
#include <vector>

#include <GridKit/AutomaticDifferentiation/Enzyme/CooEntries.hpp>
#include <GridKit/Constants.hpp>
#include <GridKit/Model/ConfigurationChecks.hpp>
#include <GridKit/ScalarTraits.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    using Log = ::GridKit::Utilities::Logger;

    /// Limit of a variable or constraint function without one
    template <typename RealT>
    inline constexpr RealT UNBOUNDED = std::numeric_limits<RealT>::infinity();

    /**
     * @brief Optimal power flow component
     *
     * Internal constraints limit the internal variables or the internal
     * constraint functions. External constraint functions add the active and
     * reactive power into each terminal bus to its balance. Local variables
     * are the internal variables followed by the real and imaginary voltage
     * of each terminal bus, and local constraint functions are the internal
     * ones followed by the external ones. The system model maps local indices
     * to global ones.
     */
    template <typename scalar_type, typename index_type>
    class Component
    {
    public:
      using ScalarT     = scalar_type;
      using IdxT        = index_type;
      using RealT       = typename ScalarTraits<ScalarT>::RealT;
      using CooEntriesT = Enzyme::Sparse::CooEntries;

      /// Local variables and constraint functions per terminal
      static constexpr IdxT TERMINAL_SIZE = 2;

      virtual ~Component() = default;

      virtual Model::ConfigurationChecks verify() const = 0;

      /**
       * @brief Set the starting point of the internal variables
       *
       * @param[out] x - Variables
       */
      virtual int initialize(ScalarT* x) = 0;

      /// Add the objective contribution at `x` to `f`
      virtual int evaluateObjective(const ScalarT* x, ScalarT& f) = 0;

      /// Add the objective gradient at `x` to `gradient`
      virtual int evaluateGradient(const ScalarT* x, RealT* gradient) = 0;

      /// Add the constraint contributions at `x` to `g`
      virtual int evaluateConstraints(const ScalarT* x, ScalarT* g) = 0;

      /// Evaluate the constraint Jacobian entries at `x` into `jacobian()`
      virtual int evaluateJacobian(const ScalarT* x) = 0;

      /// Evaluate the lower Lagrangian Hessian entries at `x` into `hessian()`
      virtual int evaluateHessian(const ScalarT* x, RealT sigma, const RealT* lambda) = 0;

      IdxT size() const
      {
        return static_cast<IdxT>(variable_indices_.size());
      }

      IdxT sizeConstraints() const
      {
        return static_cast<IdxT>(constraint_indices_.size());
      }

      IdxT sizeInternalVariables() const
      {
        return static_cast<IdxT>(variable_lower_.size());
      }

      IdxT sizeInternalConstraints() const
      {
        return static_cast<IdxT>(constraint_lower_.size());
      }

      const std::string& id() const
      {
        return id_;
      }

      bool inService() const
      {
        return in_service_ > ZERO<RealT>;
      }

      /// Bus number of each terminal
      const std::vector<IdxT>& terminals() const
      {
        return terminals_;
      }

      /// Lower limits of the internal variables
      const std::vector<RealT>& variableLower() const
      {
        return variable_lower_;
      }

      /// Upper limits of the internal variables
      const std::vector<RealT>& variableUpper() const
      {
        return variable_upper_;
      }

      /// Lower limits of the internal constraint functions
      const std::vector<RealT>& constraintLower() const
      {
        return constraint_lower_;
      }

      /// Upper limits of the internal constraint functions
      const std::vector<RealT>& constraintUpper() const
      {
        return constraint_upper_;
      }

      std::vector<IdxT>& variableIndices()
      {
        return variable_indices_;
      }

      const std::vector<IdxT>& variableIndices() const
      {
        return variable_indices_;
      }

      std::vector<IdxT>& constraintIndices()
      {
        return constraint_indices_;
      }

      const std::vector<IdxT>& constraintIndices() const
      {
        return constraint_indices_;
      }

      const CooEntriesT& jacobian() const
      {
        return jacobian_;
      }

      const CooEntriesT& hessian() const
      {
        return hessian_;
      }

    protected:
      /**
       * @param[in] id - Device `id`
       * @param[in] terminals - Bus number of each terminal
       * @param[in] in_service - Service status
       * @param[in] size_internal_variables - Number of internal variables
       * @param[in] size_internal_constraints - Number of internal constraint functions
       */
      Component(std::string       id,
                std::vector<IdxT> terminals,
                bool              in_service,
                IdxT              size_internal_variables,
                IdxT              size_internal_constraints)
        : id_(std::move(id)),
          terminals_(std::move(terminals)),
          variable_lower_(size_internal_variables, -UNBOUNDED<RealT>),
          variable_upper_(size_internal_variables, UNBOUNDED<RealT>),
          constraint_lower_(size_internal_constraints, -UNBOUNDED<RealT>),
          constraint_upper_(size_internal_constraints, UNBOUNDED<RealT>),
          variable_indices_(size_internal_variables + TERMINAL_SIZE * terminals_.size(), INVALID_INDEX<IdxT>),
          constraint_indices_(size_internal_constraints + TERMINAL_SIZE * terminals_.size(), INVALID_INDEX<IdxT>)
      {
        if (!in_service)
        {
          in_service_ = ZERO<RealT>;
        }
      }

      std::string        id_;
      std::vector<IdxT>  terminals_;
      RealT              in_service_{ONE<RealT>}; ///< \f$u\f$ Service status, 1 in service and 0 out of service
      std::vector<RealT> variable_lower_;
      std::vector<RealT> variable_upper_;
      std::vector<RealT> constraint_lower_;
      std::vector<RealT> constraint_upper_;
      std::vector<IdxT>  variable_indices_;
      std::vector<IdxT>  constraint_indices_;
      CooEntriesT        jacobian_;
      CooEntriesT        hessian_;
    };
  } // namespace OptimalPowerFlow
} // namespace GridKit
