/**
 * @file Generator.hpp
 * @brief Declaration of an optimal power flow generator.
 */

#pragma once

#include <array>
#include <cstddef>

#include <GridKit/Model/OptimalPowerFlow/Bus/Bus.hpp>
#include <GridKit/Model/OptimalPowerFlow/Component.hpp>
#include <GridKit/Model/OptimalPowerFlow/Generator/GeneratorData.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /// Internal variables of a `Generator`
    enum class GeneratorInternalVariables : size_t
    {
      P, ///< \f$P\f$ Active power into the bus [p.u.]
      Q, ///< \f$Q\f$ Reactive power into the bus [p.u.]
    };

    /// External equality contributions of a `Generator`, the power into its bus
    enum class GeneratorExternalEqualities : size_t
    {
      DP, ///< \f$\Delta P^\mathrm{bus}\f$ Active power into the bus [p.u.]
      DQ, ///< \f$\Delta Q^\mathrm{bus}\f$ Reactive power into the bus [p.u.]
    };

    /**
     * @brief Dispatchable active and reactive power injection with a
     * quadratic cost
     */
    template <typename scalar_type, typename index_type>
    class Generator final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;
      using ModelDataT = GeneratorData<RealT, IdxT>;
      using BusT       = Bus<scalar_type, index_type>;

      using ComponentT::x_lower_;
      using ComponentT::x_upper_;

    public:
      using InternalVariablesT  = GeneratorInternalVariables;
      using ExternalVariablesT  = Utilities::EmptyEnum;
      using InternalEqualitiesT = Utilities::EmptyEnum;
      using ExternalEqualitiesT = GeneratorExternalEqualities;
      using InequalitiesT       = Utilities::EmptyEnum;

      static constexpr std::array<IdxT, Utilities::enum_size<GeneratorPorts>()> NUM_PORT_VARIABLES{0};
      static constexpr std::array<IdxT, Utilities::enum_size<GeneratorPorts>()> NUM_PORT_EQUALITIES{
          Utilities::enum_size<typename BusT::InternalEqualitiesT>()};

      explicit Generator(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        initialize() override;

      int evaluateGradient() override;
      int evaluateEqualityJacobian() override;
      int evaluateInequalityJacobian() override;
      int evaluateHessian(RealT sigma, const RealT* lambda_c, const RealT* lambda_d) override;

      [[gnu::always_inline]] ScalarT objective(const ScalarT* x, const ScalarT* external) const override;
      [[gnu::always_inline]] void    externalEqualities(const ScalarT* x, const ScalarT* external, ScalarT* c) const override;

    private:
      void initializeParameters(const ModelDataT& data);
      void setDerivedParameters();

      RealT Pmin_{-UNBOUNDED<RealT>};
      RealT Pmax_{UNBOUNDED<RealT>};
      RealT Qmin_{-UNBOUNDED<RealT>};
      RealT Qmax_{UNBOUNDED<RealT>};
      RealT c0_{ZERO<RealT>};
      RealT c1_{ZERO<RealT>};
      RealT c2_{ZERO<RealT>};
      bool  in_service_{true};

      RealT u_{ONE<RealT>}; ///< Service status as a factor
    };
  } // namespace OptimalPowerFlow
} // namespace GridKit
