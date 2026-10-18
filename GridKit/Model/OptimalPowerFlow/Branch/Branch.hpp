/**
 * @file Branch.hpp
 * @brief Declaration of an optimal power flow branch.
 */

#pragma once

#include <array>

#include <GridKit/Model/OptimalPowerFlow/Branch/BranchData.hpp>
#include <GridKit/Model/OptimalPowerFlow/Bus/Bus.hpp>
#include <GridKit/Model/OptimalPowerFlow/Component.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /**
     * @brief Line or off-nominal transformer between two buses
     *
     * Power is positive into the buses.
     */
    template <typename scalar_type, typename index_type>
    class Branch final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using BusT       = Bus<scalar_type, index_type>;
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;
      using ModelDataT = BranchData<RealT, IdxT>;

      using ComponentT::d_upper_;

    public:
      static constexpr IdxT NUM_PORTS                 = 2;
      static constexpr IdxT NUM_INTERNAL_VARIABLES    = 0;
      static constexpr IdxT NUM_INTERNAL_EQUALITIES   = 0;
      static constexpr IdxT NUM_INTERNAL_INEQUALITIES = 2;

      static constexpr std::array<IdxT, NUM_PORTS> NUM_PORT_VARIABLES{BusT::NUM_INTERNAL_VARIABLES,
                                                                      BusT::NUM_INTERNAL_VARIABLES};
      static constexpr std::array<IdxT, NUM_PORTS> NUM_PORT_EQUALITIES{BusT::NUM_INTERNAL_EQUALITIES,
                                                                       BusT::NUM_INTERNAL_EQUALITIES};

      explicit Branch(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        allocate() override;

      int evaluateEqualityJacobian() override;
      int evaluateInequalityJacobian() override;
      int evaluateHessian(RealT sigma, const RealT* lambda_c, const RealT* lambda_d) override;

      [[gnu::always_inline]] void externalEqualities(const ScalarT* x, ScalarT* c) const override;
      [[gnu::always_inline]] void internalInequalities(const ScalarT* x, ScalarT* d) const override;

      bool inService() const
      {
        return in_service_ > ZERO<RealT>;
      }

    private:
      [[gnu::always_inline]] void power(const ScalarT* x, ScalarT* s) const;
      void                        setDerivedParams();

      RealT in_service_{ONE<RealT>};
      RealT R_{0.0};
      RealT X_{0.0};
      RealT G_{0.0};
      RealT B_{0.0};
      RealT Gmag_{0.0};
      RealT Bmag_{0.0};
      RealT tap_{1.0};
      RealT phase_{0.0};
      RealT Smax_{0.0};

      RealT g11_{0.0};
      RealT b11_{0.0};
      RealT g12_{0.0};
      RealT b12_{0.0};
      RealT g21_{0.0};
      RealT b21_{0.0};
      RealT g22_{0.0};
      RealT b22_{0.0};
    };
  } // namespace OptimalPowerFlow
} // namespace GridKit
