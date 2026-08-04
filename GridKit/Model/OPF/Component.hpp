#pragma once

#include <span>

#include <GridKit/Model/StateData.hpp>
#include <GridKit/Optimization/DerivativeStructure.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{
  namespace OPF
  {
    /**
     * @brief Component-local contribution to the AC optimal power flow model.
     *
     * Components own fixed local derivative structures and numerical buffers.
     * Their index spans map local variables and active local outputs into the
     * system-level nonlinear program.
     */
    template <class scalar_type, typename index_type>
    class Component
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename GridKit::ScalarTraits<ScalarT>::RealT;
      using JacobianEntryT =
          GridKit::Optimization::LocalJacobianEntry<IdxT>;
      using HessianEntryT =
          GridKit::Optimization::LocalHessianEntry<IdxT>;

      virtual ~Component() = default;

      virtual std::span<const IdxT> variableIndices() const   = 0;
      virtual std::span<const IdxT> constraintIndices() const = 0;

      virtual std::span<const JacobianEntryT> jacobianEntries() const = 0;
      virtual std::span<const HessianEntryT>  hessianEntries() const  = 0;

      virtual int gatherVariables(const ScalarT* global_variables) = 0;
      virtual int evaluateObjective()                              = 0;
      virtual int evaluateObjectiveGradient()                      = 0;
      virtual int evaluateConstraints()                            = 0;
      virtual int evaluateJacobian()                               = 0;
      virtual int evaluateHessian(
          RealT                  objective_factor,
          std::span<const RealT> local_multipliers) = 0;

      virtual RealT                    objective() const               = 0;
      virtual std::span<const ScalarT> objectiveGradientValues() const = 0;
      virtual std::span<const ScalarT> constraintValues() const        = 0;
      virtual std::span<const RealT>   jacobianValues() const          = 0;
      virtual std::span<const RealT>   hessianValues() const           = 0;

      virtual bool hasJacobian() const = 0;
      virtual bool hasHessian() const  = 0;

      virtual int initialize([[maybe_unused]] ScalarT* global_variables) const
      {
        return 0;
      }

      virtual int setVariableBounds(
          [[maybe_unused]] RealT* global_lower_bounds,
          [[maybe_unused]] RealT* global_upper_bounds) const
      {
        return 0;
      }

      virtual int updateSolutionState(
          [[maybe_unused]] const ScalarT*    global_variables,
          [[maybe_unused]] Model::StateData& state) const
      {
        return 0;
      }
    };
  } // namespace OPF
} // namespace GridKit
