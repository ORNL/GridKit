#pragma once

#include <array>
#include <cstddef>
#include <type_traits>

#include <GridKit/Model/OptimalPowerFlow/ComponentImpl.hpp>

namespace GridKit
{
  namespace Testing
  {
    namespace OPF = GridKit::OptimalPowerFlow;

    /// Two variables and two equality constraints for other components to read and add to
    class EqualityOwner final : public OPF::Component<double, size_t>
    {
      using RealT = double;
      using IdxT  = size_t;

    public:
      enum class Variables : size_t
      {
        A,
        B,
      };

      enum class Constraints : size_t
      {
        A,
        B,
      };

      using InternalVariablesT   = Variables;
      using ExternalVariablesT   = Utilities::EmptyEnum;
      using InternalConstraintsT = Constraints;
      using ExternalConstraintsT = Utilities::EmptyEnum;

      static constexpr std::array<IdxT, 0> NUM_PORT_VARIABLES{};
      static constexpr std::array<IdxT, 0> NUM_PORT_CONSTRAINTS{};

      EqualityOwner()
        : Component(std::type_identity<EqualityOwner>{}, "owner")
      {
      }

      Model::ConfigurationChecks verify() const override
      {
        return {};
      }

      void internalConstraints(const RealT*, const RealT*, RealT* g) const override
      {
        g[0] = 10.0;
        g[1] = -10.0;
      }
    };

    /// One internal variable, equality, and inequality, coupled to an owner at three ports
    class CoupledComponent final : public OPF::Component<double, size_t>
    {
      using RealT = double;
      using IdxT  = size_t;

    public:
      enum class InternalVariables : size_t
      {
        X,
      };

      enum class ExternalVariables : size_t
      {
        Y, ///< Read at port 0
        Z, ///< Read at port 1
      };

      enum class InternalConstraints : size_t
      {
        C, ///< Equality
        D, ///< Inequality
      };

      enum class ExternalConstraints : size_t
      {
        C0, ///< Added at port 0
        C2, ///< Added at port 2
      };

      using InternalVariablesT   = InternalVariables;
      using ExternalVariablesT   = ExternalVariables;
      using InternalConstraintsT = InternalConstraints;
      using ExternalConstraintsT = ExternalConstraints;

      static constexpr std::array<IdxT, 3> NUM_PORT_VARIABLES{1, 1, 0};
      static constexpr std::array<IdxT, 3> NUM_PORT_CONSTRAINTS{1, 0, 1};

      CoupledComponent()
        : Component(std::type_identity<CoupledComponent>{}, "coupled")
      {
      }

      Model::ConfigurationChecks verify() const override
      {
        return {};
      }

      int tagInequality() override
      {
        tag_[static_cast<size_t>(InternalConstraints::C)] = false;
        tag_[static_cast<size_t>(InternalConstraints::D)] = true;
        return 0;
      }

      int evaluateGradient() override;
      int evaluateJacobian() override;
      int evaluateHessian(RealT sigma) override;

      [[gnu::always_inline]] RealT objective(const RealT* x, const RealT* external) const override
      {
        return 0.5 * x[0] * x[0] + external[0] * external[1];
      }

      [[gnu::always_inline]] void internalConstraints(const RealT* x, const RealT* external, RealT* g) const override
      {
        g[0] = x[0] * external[0];
        g[1] = x[0] * x[0] + external[0] * external[1];
      }

      [[gnu::always_inline]] void externalConstraints(const RealT* x, const RealT* external, RealT* g) const override
      {
        g[0] = external[0] * external[1];
        g[1] = x[0] * external[1];
      }
    };
  } // namespace Testing
} // namespace GridKit
