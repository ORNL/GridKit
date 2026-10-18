/**
 * @file Generator.hpp
 * @brief Declaration of an optimal power flow generator.
 */

#pragma once

#include <array>

#include <GridKit/Model/OptimalPowerFlow/Bus/Bus.hpp>
#include <GridKit/Model/OptimalPowerFlow/Component.hpp>
#include <GridKit/Model/OptimalPowerFlow/Generator/GeneratorData.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /**
     * @brief Dispatchable active and reactive power injection with a
     * quadratic cost
     */
    template <typename scalar_type, typename index_type>
    class Generator final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using BusT       = Bus<scalar_type, index_type>;
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;
      using ModelDataT = GeneratorData<RealT, IdxT>;

      using ComponentT::x_;
      using ComponentT::x_lower_;
      using ComponentT::x_upper_;

    public:
      static constexpr IdxT NUM_PORTS                 = 1;
      static constexpr IdxT NUM_INTERNAL_VARIABLES    = 2;
      static constexpr IdxT NUM_INTERNAL_EQUALITIES   = 0;
      static constexpr IdxT NUM_INTERNAL_INEQUALITIES = 0;

      static constexpr std::array<IdxT, NUM_PORTS> NUM_PORT_VARIABLES{0};
      static constexpr std::array<IdxT, NUM_PORTS> NUM_PORT_EQUALITIES{BusT::NUM_INTERNAL_EQUALITIES};

      explicit Generator(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        allocate() override;
      int                        initialize() override;

      int evaluateGradient() override;
      int evaluateEqualityJacobian() override;
      int evaluateHessian(RealT sigma, const RealT* lambda_c, const RealT* lambda_d) override;

      [[gnu::always_inline]] ScalarT objective(const ScalarT* x) const override;
      [[gnu::always_inline]] void    externalEqualities(const ScalarT* x, ScalarT* c) const override;

    private:
      RealT in_service_{ONE<RealT>};
      RealT Pmin_{0.0};
      RealT Pmax_{0.0};
      RealT Qmin_{0.0};
      RealT Qmax_{0.0};
      RealT c0_{0.0};
      RealT c1_{0.0};
      RealT c2_{0.0};
    };
  } // namespace OptimalPowerFlow
} // namespace GridKit
