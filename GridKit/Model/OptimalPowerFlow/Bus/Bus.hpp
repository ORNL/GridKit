/**
 * @file Bus.hpp
 * @brief Declaration of an optimal power flow bus.
 */

#pragma once

#include <array>

#include <GridKit/Model/OptimalPowerFlow/Bus/BusData.hpp>
#include <GridKit/Model/OptimalPowerFlow/Component.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /**
     * @brief Bus voltage in Cartesian coordinates, power balance, and shunt
     * admittance to ground
     *
     * A port to the bus reads its voltage, the internal variables, and
     * contributes to its power balance, the internal equality rows. Each
     * balance row is the power into the bus from its shunt and from its ports.
     */
    template <typename scalar_type, typename index_type>
    class Bus final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;
      using ModelDataT = BusData<RealT, IdxT>;

      using ComponentT::d_lower_;
      using ComponentT::d_upper_;
      using ComponentT::x_;
      using ComponentT::x_lower_;
      using ComponentT::x_upper_;

    public:
      static constexpr IdxT NUM_PORTS                 = 0;
      static constexpr IdxT NUM_INTERNAL_VARIABLES    = 2;
      static constexpr IdxT NUM_INTERNAL_EQUALITIES   = 2;
      static constexpr IdxT NUM_INTERNAL_INEQUALITIES = 1;

      static constexpr std::array<IdxT, NUM_PORTS> NUM_PORT_VARIABLES{};
      static constexpr std::array<IdxT, NUM_PORTS> NUM_PORT_EQUALITIES{};

      explicit Bus(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        allocate() override;
      int                        initialize() override;

      int evaluateEqualityJacobian() override;
      int evaluateInequalityJacobian() override;
      int evaluateHessian(RealT sigma, const RealT* lambda_c, const RealT* lambda_d) override;

      [[gnu::always_inline]] void internalEqualities(const ScalarT* x, ScalarT* c) const override;
      [[gnu::always_inline]] void internalInequalities(const ScalarT* x, ScalarT* d) const override;

      /// Angle reference, with zero imaginary voltage
      bool reference() const
      {
        return reference_;
      }

    private:
      bool  reference_{false};
      RealT Vmin_{0.0};
      RealT Vmax_{0.0};
      RealT G_{0.0};
      RealT B_{0.0};
    };
  } // namespace OptimalPowerFlow
} // namespace GridKit
