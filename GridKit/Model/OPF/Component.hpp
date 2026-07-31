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
     * System owns global storage. Components own semantic variable bindings
     * and add their objective, constraint, and derivative contributions to
     * that storage.
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
      virtual IdxT nnz() const                     = 0;

      virtual void setVariableOffset([[maybe_unused]] IdxT offset)
      {
        if (sizeInternalVariables() != 0)
        {
          throw std::logic_error("OPF component does not bind its internal variables");
        }
      }

      virtual void setConstraintOffset([[maybe_unused]] IdxT offset)
      {
        if (sizeInternalConstraints() != 0)
        {
          throw std::logic_error("OPF component does not bind its internal constraints");
        }
      }

      virtual void addJacobianSparsity(
          [[maybe_unused]] std::vector<JacobianEntry>& entries) const
      {
        if (nnz() != 0)
        {
          throw std::logic_error("OPF component does not define Jacobian sparsity");
        }
      }

      virtual int setJacobianSlots(const std::vector<IdxT>& slots)
      {
        if (slots.empty() && nnz() == 0)
        {
          return 0;
        }
        return 1;
      }

      virtual int initialize([[maybe_unused]] ScalarT* variables) const
      {
        return 0;
      }

      virtual int setVariableBounds([[maybe_unused]] RealT* lower,
                                    [[maybe_unused]] RealT* upper) const
      {
        return 0;
      }

      virtual int setConstraintBounds([[maybe_unused]] RealT* lower,
                                      [[maybe_unused]] RealT* upper) const
      {
        return 0;
      }

      virtual int addObjective([[maybe_unused]] const ScalarT* variables,
                               [[maybe_unused]] RealT&         value) const
      {
        return 0;
      }

      virtual int addObjectiveGradient([[maybe_unused]] const ScalarT* variables,
                                       [[maybe_unused]] ScalarT*       gradient) const
      {
        return 0;
      }

      virtual int addConstraints([[maybe_unused]] const ScalarT* variables,
                                 [[maybe_unused]] ScalarT*       constraints) const
      {
        return 0;
      }

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
