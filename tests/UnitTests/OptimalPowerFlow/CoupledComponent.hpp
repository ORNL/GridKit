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

    /// Two variables and two equality rows for other components to read and add to
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

      enum class Equalities : size_t
      {
        A,
        B,
      };

      using InternalVariablesT  = Variables;
      using ExternalVariablesT  = Utilities::EmptyEnum;
      using InternalEqualitiesT = Equalities;
      using ExternalEqualitiesT = Utilities::EmptyEnum;
      using InequalitiesT       = Utilities::EmptyEnum;

      static constexpr std::array<IdxT, 0> NUM_PORT_VARIABLES{};
      static constexpr std::array<IdxT, 0> NUM_PORT_EQUALITIES{};

      EqualityOwner()
        : Component(std::type_identity<EqualityOwner>{}, "owner")
      {
      }

      Model::ConfigurationChecks verify() const override
      {
        return {};
      }

      void internalEqualities(const RealT*, const RealT*, RealT* c) const override
      {
        c[0] = 10.0;
        c[1] = -10.0;
      }
    };

    /// One internal variable, row, and inequality, coupled to an owner at three ports
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

      enum class InternalEqualities : size_t
      {
        C,
      };

      enum class ExternalEqualities : size_t
      {
        C0, ///< Added at port 0
        C2, ///< Added at port 2
      };

      enum class Inequalities : size_t
      {
        D,
      };

      using InternalVariablesT  = InternalVariables;
      using ExternalVariablesT  = ExternalVariables;
      using InternalEqualitiesT = InternalEqualities;
      using ExternalEqualitiesT = ExternalEqualities;
      using InequalitiesT       = Inequalities;

      static constexpr std::array<IdxT, 3> NUM_PORT_VARIABLES{1, 1, 0};
      static constexpr std::array<IdxT, 3> NUM_PORT_EQUALITIES{1, 0, 1};

      CoupledComponent()
        : Component(std::type_identity<CoupledComponent>{}, "coupled")
      {
      }

      Model::ConfigurationChecks verify() const override
      {
        return {};
      }

      int evaluateGradient() override;
      int evaluateEqualityJacobian() override;
      int evaluateInequalityJacobian() override;
      int evaluateHessian(RealT sigma, const RealT* lambda_c, const RealT* lambda_d) override;

      [[gnu::always_inline]] RealT objective(const RealT* x, const RealT* external) const override
      {
        return 0.5 * x[0] * x[0] + external[0] * external[1];
      }

      [[gnu::always_inline]] void internalEqualities(const RealT* x, const RealT* external, RealT* c) const override
      {
        c[0] = x[0] * external[0];
      }

      [[gnu::always_inline]] void externalEqualities(const RealT* x, const RealT* external, RealT* c) const override
      {
        c[0] = external[0] * external[1];
        c[1] = x[0] * external[1];
      }

      [[gnu::always_inline]] void inequalities(const RealT* x, const RealT* external, RealT* d) const override
      {
        d[0] = x[0] * x[0] + external[0] * external[1];
      }
    };
  } // namespace Testing
} // namespace GridKit
