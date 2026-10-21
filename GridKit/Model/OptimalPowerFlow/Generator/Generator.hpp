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

    /// External constraint contributions of a `Generator`, the power into its bus
    enum class GeneratorExternalConstraints : size_t
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

      using ComponentT::x_;
      using ComponentT::x_lower_;
      using ComponentT::x_upper_;

    public:
      using InternalVariablesT   = GeneratorInternalVariables;
      using ExternalVariablesT   = Utilities::EmptyEnum;
      using InternalConstraintsT = Utilities::EmptyEnum;
      using ExternalConstraintsT = GeneratorExternalConstraints;

      static constexpr std::array<IdxT, Utilities::enum_size<GeneratorPorts>()> NUM_PORT_VARIABLES{0};
      static constexpr std::array<IdxT, Utilities::enum_size<GeneratorPorts>()> NUM_PORT_CONSTRAINTS{
          BusT::NUM_SHARED_CONSTRAINTS};

      explicit Generator(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        initialize() override;

      int evaluateGradient() override;
      int evaluateJacobian() override;
      int evaluateHessian(RealT sigma) override;

      [[gnu::always_inline]] ScalarT objective(const ScalarT* x, const ScalarT* external) const override;
      [[gnu::always_inline]] void    externalConstraints(const ScalarT* x, const ScalarT* external, ScalarT* g) const override;

      /// Active power into the bus [p.u.]
      const ScalarT& P() const
      {
        return x_.getData()[static_cast<size_t>(GeneratorInternalVariables::P)];
      }

      /// Reactive power into the bus [p.u.]
      const ScalarT& Q() const
      {
        return x_.getData()[static_cast<size_t>(GeneratorInternalVariables::Q)];
      }

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
