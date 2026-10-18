/**
 * @file Load.hpp
 * @brief Declaration of an optimal power flow load.
 */

#pragma once

#include <array>
#include <cstddef>

#include <GridKit/Model/OptimalPowerFlow/Bus/Bus.hpp>
#include <GridKit/Model/OptimalPowerFlow/Component.hpp>
#include <GridKit/Model/OptimalPowerFlow/Load/LoadData.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /// External variables of a `Load`, the voltage of its bus
    enum class LoadExternalVariables : size_t
    {
      VR, ///< \f$V_\mathrm{r}\f$ Real voltage [p.u.]
      VI, ///< \f$V_\mathrm{i}\f$ Imaginary voltage [p.u.]
    };

    /// External equality contributions of a `Load`, the power into its bus
    enum class LoadExternalEqualities : size_t
    {
      DP, ///< \f$\Delta P^\mathrm{bus}\f$ Active power into the bus [p.u.]
      DQ, ///< \f$\Delta Q^\mathrm{bus}\f$ Reactive power into the bus [p.u.]
    };

    /**
     * @brief ZIP demand: constant-power, constant-current, and
     * constant-impedance parts
     */
    template <typename scalar_type, typename index_type>
    class Load final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;
      using ModelDataT = LoadData<RealT, IdxT>;
      using BusT       = Bus<scalar_type, index_type>;

    public:
      using InternalVariablesT  = Utilities::EmptyEnum;
      using ExternalVariablesT  = LoadExternalVariables;
      using InternalEqualitiesT = Utilities::EmptyEnum;
      using ExternalEqualitiesT = LoadExternalEqualities;
      using InequalitiesT       = Utilities::EmptyEnum;

      static constexpr std::array<IdxT, Utilities::enum_size<LoadPorts>()> NUM_PORT_VARIABLES{
          Utilities::enum_size<typename BusT::InternalVariablesT>()};
      static constexpr std::array<IdxT, Utilities::enum_size<LoadPorts>()> NUM_PORT_EQUALITIES{
          Utilities::enum_size<typename BusT::InternalEqualitiesT>()};

      explicit Load(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;

      int evaluateGradient() override;
      int evaluateEqualityJacobian() override;
      int evaluateInequalityJacobian() override;
      int evaluateHessian(RealT sigma, const RealT* lambda_c, const RealT* lambda_d) override;

      [[gnu::always_inline]] void externalEqualities(const ScalarT* x, const ScalarT* external, ScalarT* c) const override;

    private:
      void initializeParameters(const ModelDataT& data);
      void setDerivedParameters();

      RealT P_{ZERO<RealT>};
      RealT Q_{ZERO<RealT>};
      RealT Pi_{ZERO<RealT>};
      RealT Qi_{ZERO<RealT>};
      RealT Pz_{ZERO<RealT>};
      RealT Qz_{ZERO<RealT>};
      bool  in_service_{true};

      RealT u_{ONE<RealT>};           ///< Service status as a factor
      RealT current_on_{ZERO<RealT>}; ///< 1 with an in-service constant-current part, which needs the voltage magnitude
    };
  } // namespace OptimalPowerFlow
} // namespace GridKit
