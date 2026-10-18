/**
 * @file Load.hpp
 * @brief Declaration of an optimal power flow load.
 */

#pragma once

#include <array>

#include <GridKit/Model/OptimalPowerFlow/Bus/Bus.hpp>
#include <GridKit/Model/OptimalPowerFlow/Component.hpp>
#include <GridKit/Model/OptimalPowerFlow/Load/LoadData.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /**
     * @brief ZIP demand: constant-power, constant-current, and
     * constant-impedance parts
     */
    template <typename scalar_type, typename index_type>
    class Load final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using BusT       = Bus<scalar_type, index_type>;
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;
      using ModelDataT = LoadData<RealT, IdxT>;

    public:
      static constexpr IdxT NUM_PORTS                 = 1;
      static constexpr IdxT NUM_INTERNAL_VARIABLES    = 0;
      static constexpr IdxT NUM_INTERNAL_EQUALITIES   = 0;
      static constexpr IdxT NUM_INTERNAL_INEQUALITIES = 0;

      static constexpr std::array<IdxT, NUM_PORTS> NUM_PORT_VARIABLES{BusT::NUM_INTERNAL_VARIABLES};
      static constexpr std::array<IdxT, NUM_PORTS> NUM_PORT_EQUALITIES{BusT::NUM_INTERNAL_EQUALITIES};

      explicit Load(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;

      int evaluateEqualityJacobian() override;
      int evaluateHessian(RealT sigma, const RealT* lambda_c, const RealT* lambda_d) override;

      [[gnu::always_inline]] void externalEqualities(const ScalarT* x, ScalarT* c) const override;

    private:
      RealT in_service_{ONE<RealT>};
      RealT current_part_{ZERO<RealT>}; ///< 1 with a constant-current part, which needs the voltage magnitude, else 0
      RealT P_{0.0};
      RealT Q_{0.0};
      RealT Pi_{0.0};
      RealT Qi_{0.0};
      RealT Pz_{0.0};
      RealT Qz_{0.0};
    };
  } // namespace OptimalPowerFlow
} // namespace GridKit
