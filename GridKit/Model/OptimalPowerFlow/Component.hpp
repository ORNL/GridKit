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
#include <GridKit/Model/StateData.hpp>
#include <GridKit/ScalarTraits.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    using Log = ::GridKit::Utilities::Logger;

    /// Bound of a variable or constraint without a limit
    template <typename RealT>
    inline constexpr RealT UNBOUNDED = std::numeric_limits<RealT>::infinity();

    /**
     * @brief Optimal power flow component
     *
     * Local variables are the internal variables followed by the real and
     * imaginary voltage of each terminal bus. Local constraints are the
     * internal constraints followed by the active and reactive power into
     * each terminal bus. The system model maps local indices to global ones.
     */
    template <typename scalar_type, typename index_type>
    class Component
    {
    public:
      using ScalarT     = scalar_type;
      using IdxT        = index_type;
      using RealT       = typename ScalarTraits<ScalarT>::RealT;
      using CooEntriesT = Enzyme::Sparse::CooEntries;

      /// Local variables and constraints per terminal
      static constexpr IdxT TERMINAL_SIZE = 2;

      virtual ~Component() = default;

      virtual int verify() const = 0;

      /**
       * @brief Set initial values and bounds of the internal variables
       *
       * @param[in] state - Operating state
       * @param[out] x - Variables
       * @param[out] x_lower - Variable lower bounds
       * @param[out] x_upper - Variable upper bounds
       */
      virtual int initialize(const Model::StateData& state,
                             ScalarT*                x,
                             RealT*                  x_lower,
                             RealT*                  x_upper) = 0;

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

      /// Active and reactive power into each terminal bus at `x`
      virtual int evaluateTerminalPower(const ScalarT* x, ScalarT* power) = 0;

      IdxT size() const
      {
        return static_cast<IdxT>(variable_indices_.size());
      }

      IdxT sizeConstraints() const
      {
        return static_cast<IdxT>(constraint_indices_.size());
      }

      IdxT sizeInternal() const
      {
        return size_internal_;
      }

      IdxT sizeInternalConstraints() const
      {
        return static_cast<IdxT>(constraint_lower_.size());
      }

      const std::string& id() const
      {
        return id_;
      }

      /// Bus number of each terminal
      const std::vector<IdxT>& terminals() const
      {
        return terminals_;
      }

      /// Lower bounds of the internal constraints
      const std::vector<RealT>& constraintLower() const
      {
        return constraint_lower_;
      }

      /// Upper bounds of the internal constraints
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
       * @param[in] id - Case device `id`
       * @param[in] terminals - Bus number of each terminal
       * @param[in] size_internal - Number of internal variables
       * @param[in] size_internal_constraints - Number of internal constraints
       */
      Component(std::string       id,
                std::vector<IdxT> terminals,
                IdxT              size_internal,
                IdxT              size_internal_constraints)
        : id_(std::move(id)),
          terminals_(std::move(terminals)),
          size_internal_(size_internal),
          constraint_lower_(size_internal_constraints, -UNBOUNDED<RealT>),
          constraint_upper_(size_internal_constraints, UNBOUNDED<RealT>),
          variable_indices_(size_internal + TERMINAL_SIZE * terminals_.size(), INVALID_INDEX<IdxT>),
          constraint_indices_(size_internal_constraints + TERMINAL_SIZE * terminals_.size(), INVALID_INDEX<IdxT>)
      {
      }

      /**
       * @brief Log `message` for this component if `condition` fails
       *
       * @return 0 if `condition` holds, 1 otherwise
       */
      int check(bool condition, const char* message) const
      {
        if (condition)
        {
          return 0;
        }
        Log::error() << id_ << ": " << message << "\n";
        return 1;
      }

      /**
       * @brief Power into the bus at `terminal` given by the state
       *
       * @return false if the state has no current for the terminal or no
       * voltage for its bus
       */
      bool statePower(const Model::StateData& state, IdxT terminal, RealT& p, RealT& q) const
      {
        return Model::terminalPower(state, id_, terminals_[terminal], terminal, terminals_.size(), p, q);
      }

      std::string        id_;
      std::vector<IdxT>  terminals_;
      IdxT               size_internal_{0};
      std::vector<RealT> constraint_lower_;
      std::vector<RealT> constraint_upper_;
      std::vector<IdxT>  variable_indices_;
      std::vector<IdxT>  constraint_indices_;
      CooEntriesT        jacobian_;
      CooEntriesT        hessian_;
    };
  } // namespace OptimalPowerFlow
} // namespace GridKit
