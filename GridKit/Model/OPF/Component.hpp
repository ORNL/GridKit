#pragma once

#include <stdexcept>
#include <utility>
#include <vector>

#include <GridKit/Model/StateData.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{
  namespace OPF
  {
    /**
     * @brief Local contribution to the system-level nonlinear program.
     *
     * System owns global storage and uses this interface to assemble every
     * concrete component into one nonlinear program. Components own their
     * internal index ranges, bind the external indices they consume, and add
     * their objective, constraint, and derivative contributions to System's
     * storage.
     */
    template <class scalar_type, typename index_type>
    class Component
    {
    public:
      using ScalarT       = scalar_type;
      using IdxT          = index_type;
      using RealT         = typename GridKit::ScalarTraits<ScalarT>::RealT;
      using JacobianEntry = std::pair<IdxT, IdxT>;

      virtual ~Component() = default;

      virtual IdxT sizeInternalVariables() const   = 0;
      virtual IdxT sizeInternalConstraints() const = 0;

      /**
       * @brief Number of ordered local Jacobian contributions.
       *
       * This count is not necessarily the number of unique entries contributed
       * to the system CSR matrix. Contributions from multiple components may
       * share the same global row and column.
       */
      virtual IdxT nnz() const = 0;

      /// Assign the first global index of this component's internal variables.
      virtual void setVariableOffset([[maybe_unused]] IdxT offset)
      {
        if (sizeInternalVariables() != 0)
        {
          throw std::logic_error("OPF component does not bind its internal variables");
        }
      }

      /// Assign the first global row of this component's internal constraints.
      virtual void setConstraintOffset([[maybe_unused]] IdxT offset)
      {
        if (sizeInternalConstraints() != 0)
        {
          throw std::logic_error("OPF component does not bind its internal constraints");
        }
      }

      /**
       * @brief Append this component's ordered Jacobian coordinates.
       *
       * A component with a nonzero nnz() appends exactly nnz() entries. Their
       * order defines the order expected by setJacobianSlots().
       */
      virtual void addJacobianSparsity(
          [[maybe_unused]] std::vector<JacobianEntry>& entries) const
      {
        if (nnz() != 0)
        {
          throw std::logic_error("OPF component does not define Jacobian sparsity");
        }
      }

      /**
       * @brief Bind local Jacobian contributions to global CSR value indices.
       *
       * slots has one entry for every coordinate appended by
       * addJacobianSparsity(). Repeated indices are valid when multiple local
       * contributions accumulate into the same global matrix entry.
       */
      virtual int setJacobianSlots(const std::vector<IdxT>& slots)
      {
        if (slots.empty() && nnz() == 0)
        {
          return 0;
        }
        return 1;
      }

      /// Initialize only the internal variables owned by this component.
      virtual int initialize([[maybe_unused]] ScalarT* variables) const
      {
        return 0;
      }

      /// Write bounds only for the internal variables owned by this component.
      virtual int setVariableBounds([[maybe_unused]] RealT* lower,
                                    [[maybe_unused]] RealT* upper) const
      {
        return 0;
      }

      /// Write bounds only for the internal constraints owned by this component.
      virtual int setConstraintBounds([[maybe_unused]] RealT* lower,
                                      [[maybe_unused]] RealT* upper) const
      {
        return 0;
      }

      /// Accumulate this component's contribution into value.
      virtual int addObjective([[maybe_unused]] const ScalarT* variables,
                               [[maybe_unused]] RealT&         value) const
      {
        return 0;
      }

      /// Accumulate this component's objective derivatives into gradient.
      virtual int addObjectiveGradient([[maybe_unused]] const ScalarT* variables,
                                       [[maybe_unused]] ScalarT*       gradient) const
      {
        return 0;
      }

      /// Accumulate this component's equations into constraints.
      virtual int addConstraints([[maybe_unused]] const ScalarT* variables,
                                 [[maybe_unused]] ScalarT*       constraints) const
      {
        return 0;
      }

      /// Accumulate this component's derivatives into the bound CSR values.
      virtual int addJacobian([[maybe_unused]] const ScalarT* variables,
                              [[maybe_unused]] RealT*         jacobian_values) const
      {
        return 0;
      }

      /**
       * @brief Copy this component's accepted solution into portable state.
       *
       * Components without solved or voltage-dependent state preserve the
       * corresponding input-state entries by default.
       */
      virtual int updateSolutionState([[maybe_unused]] const ScalarT*    variables,
                                      [[maybe_unused]] Model::StateData& state) const
      {
        return 0;
      }
    };

  } // namespace OPF
} // namespace GridKit
