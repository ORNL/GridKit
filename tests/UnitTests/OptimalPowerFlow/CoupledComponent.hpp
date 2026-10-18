#pragma once

#include <array>
#include <type_traits>

#include <GridKit/Model/OptimalPowerFlow/ComponentImpl.hpp>

namespace GridKit
{
  namespace Testing
  {
    namespace OPF = GridKit::OptimalPowerFlow;

    class EqualityOwner final : public OPF::Component<double, size_t>
    {
      using RealT = double;
      using IdxT  = size_t;

    public:
      static constexpr IdxT                        NUM_PORTS                 = 0;
      static constexpr IdxT                        NUM_INTERNAL_VARIABLES    = 2;
      static constexpr IdxT                        NUM_INTERNAL_EQUALITIES   = 2;
      static constexpr IdxT                        NUM_INTERNAL_INEQUALITIES = 0;
      static constexpr std::array<IdxT, NUM_PORTS> NUM_PORT_VARIABLES{};
      static constexpr std::array<IdxT, NUM_PORTS> NUM_PORT_EQUALITIES{};

      EqualityOwner()
        : Component(std::type_identity<EqualityOwner>{}, "owner")
      {
      }

      Model::ConfigurationChecks verify() const override
      {
        return {};
      }

      void internalEqualities(const RealT*, RealT* c) const override
      {
        c[0] = 10.0;
        c[1] = -10.0;
      }
    };

    class CoupledComponent final : public OPF::Component<double, size_t>
    {
      using RealT = double;
      using IdxT  = size_t;

    public:
      static constexpr IdxT                        NUM_PORTS                 = 3;
      static constexpr IdxT                        NUM_INTERNAL_VARIABLES    = 1;
      static constexpr IdxT                        NUM_INTERNAL_EQUALITIES   = 1;
      static constexpr IdxT                        NUM_INTERNAL_INEQUALITIES = 1;
      static constexpr std::array<IdxT, NUM_PORTS> NUM_PORT_VARIABLES{1, 1, 0};
      static constexpr std::array<IdxT, NUM_PORTS> NUM_PORT_EQUALITIES{1, 0, 1};

      CoupledComponent()
        : Component(std::type_identity<CoupledComponent>{}, "coupled")
      {
      }

      Model::ConfigurationChecks verify() const override
      {
        return {};
      }

      [[gnu::always_inline]] RealT objective(const RealT* x) const override
      {
        return 0.5 * x[0] * x[0] + x[1] * x[2];
      }

      [[gnu::always_inline]] void internalEqualities(const RealT* x, RealT* c) const override
      {
        c[0] = x[0] * x[1];
      }

      [[gnu::always_inline]] void externalEqualities(const RealT* x, RealT* c) const override
      {
        c[0] = x[1] * x[2];
        c[1] = x[0] * x[2];
      }

      [[gnu::always_inline]] void internalInequalities(const RealT* x, RealT* d) const override
      {
        d[0] = x[0] * x[0] + x[1] * x[2];
      }

      int evaluateGradient() override;

      int evaluateEqualityJacobian() override;

      int evaluateInequalityJacobian() override;

      int evaluateHessian(RealT sigma, const RealT* lambda_c, const RealT* lambda_d) override;
    };

  } // namespace Testing
} // namespace GridKit
