/**
 * @file Branch.hpp
 * @brief Declaration of an optimal power flow branch.
 */

#pragma once

#include <array>
#include <cstddef>

#include <GridKit/Model/OptimalPowerFlow/Branch/BranchData.hpp>
#include <GridKit/Model/OptimalPowerFlow/Bus/Bus.hpp>
#include <GridKit/Model/OptimalPowerFlow/Component.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /// External variables of a `Branch`, the voltages at its ports
    enum class BranchExternalVariables : size_t
    {
      VR1, ///< \f$V_{\mathrm{r}1}\f$ Bus-1 real voltage [p.u.]
      VI1, ///< \f$V_{\mathrm{i}1}\f$ Bus-1 imaginary voltage [p.u.]
      VR2, ///< \f$V_{\mathrm{r}2}\f$ Bus-2 real voltage [p.u.]
      VI2, ///< \f$V_{\mathrm{i}2}\f$ Bus-2 imaginary voltage [p.u.]
    };

    /// External equality contributions of a `Branch`, the power into the buses at its ports
    enum class BranchExternalEqualities : size_t
    {
      DP1, ///< \f$\Delta P_1^\mathrm{bus}\f$ Active power into bus 1 [p.u.]
      DQ1, ///< \f$\Delta Q_1^\mathrm{bus}\f$ Reactive power into bus 1 [p.u.]
      DP2, ///< \f$\Delta P_2^\mathrm{bus}\f$ Active power into bus 2 [p.u.]
      DQ2, ///< \f$\Delta Q_2^\mathrm{bus}\f$ Reactive power into bus 2 [p.u.]
    };

    /// Inequality rows of a `Branch`
    enum class BranchInequalities : size_t
    {
      S1, ///< \f$P_1^2 + Q_1^2\f$ Squared apparent power at bus 1 [p.u.]
      S2, ///< \f$P_2^2 + Q_2^2\f$ Squared apparent power at bus 2 [p.u.]
    };

    /**
     * @brief Line or off-nominal transformer between two buses
     *
     * Power is positive into the buses.
     */
    template <typename scalar_type, typename index_type>
    class Branch final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;
      using ModelDataT = BranchData<RealT, IdxT>;
      using BusT       = Bus<scalar_type, index_type>;

      using ComponentT::d_upper_;

    public:
      using InternalVariablesT  = Utilities::EmptyEnum;
      using ExternalVariablesT  = BranchExternalVariables;
      using InternalEqualitiesT = Utilities::EmptyEnum;
      using ExternalEqualitiesT = BranchExternalEqualities;
      using InequalitiesT       = BranchInequalities;

      static constexpr std::array<IdxT, Utilities::enum_size<BranchPorts>()> NUM_PORT_VARIABLES{
          Utilities::enum_size<typename BusT::InternalVariablesT>(),
          Utilities::enum_size<typename BusT::InternalVariablesT>()};
      static constexpr std::array<IdxT, Utilities::enum_size<BranchPorts>()> NUM_PORT_EQUALITIES{
          Utilities::enum_size<typename BusT::InternalEqualitiesT>(),
          Utilities::enum_size<typename BusT::InternalEqualitiesT>()};

      explicit Branch(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        initialize() override;

      int evaluateGradient() override;
      int evaluateEqualityJacobian() override;
      int evaluateInequalityJacobian() override;
      int evaluateHessian(RealT sigma, const RealT* lambda_c, const RealT* lambda_d) override;

      [[gnu::always_inline]] void externalEqualities(const ScalarT* x, const ScalarT* external, ScalarT* c) const override;
      [[gnu::always_inline]] void inequalities(const ScalarT* x, const ScalarT* external, ScalarT* d) const override;

      bool inService() const
      {
        return in_service_;
      }

    private:
      void initializeParameters(const ModelDataT& data);
      void setDerivedParameters();

      [[gnu::always_inline]] void power(const ScalarT* external, ScalarT* s) const;

      RealT R_{ZERO<RealT>};
      RealT X_{ZERO<RealT>};
      RealT G_{ZERO<RealT>};
      RealT B_{ZERO<RealT>};
      RealT Gmag_{ZERO<RealT>};
      RealT Bmag_{ZERO<RealT>};
      RealT tap_{ONE<RealT>};
      RealT phase_{ZERO<RealT>};
      RealT Smax_{UNBOUNDED<RealT>};
      bool  in_service_{true};

      RealT g11_{ZERO<RealT>};
      RealT b11_{ZERO<RealT>};
      RealT g12_{ZERO<RealT>};
      RealT b12_{ZERO<RealT>};
      RealT g21_{ZERO<RealT>};
      RealT b21_{ZERO<RealT>};
      RealT g22_{ZERO<RealT>};
      RealT b22_{ZERO<RealT>};
    };
  } // namespace OptimalPowerFlow
} // namespace GridKit
