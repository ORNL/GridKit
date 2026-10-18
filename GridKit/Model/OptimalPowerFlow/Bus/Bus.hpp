/**
 * @file Bus.hpp
 * @brief Declaration of an optimal power flow bus.
 */

#pragma once

#include <array>
#include <cstddef>

#include <GridKit/Model/OptimalPowerFlow/Bus/BusData.hpp>
#include <GridKit/Model/OptimalPowerFlow/Component.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /// Internal variables of a `Bus`
    enum class BusInternalVariables : size_t
    {
      VR, ///< \f$V_\mathrm{r}\f$ Real voltage [p.u.]
      VI, ///< \f$V_\mathrm{i}\f$ Imaginary voltage [p.u.]
    };

    /// Internal equality rows of a `Bus`, its power balance
    enum class BusInternalEqualities : size_t
    {
      P, ///< Active power balance [p.u.]
      Q, ///< Reactive power balance [p.u.]
    };

    /// Inequality rows of a `Bus`
    enum class BusInequalities : size_t
    {
      VM2, ///< \f$V_\mathrm{r}^2 + V_\mathrm{i}^2\f$ Squared voltage magnitude [p.u.]
    };

    /**
     * @brief Bus voltage in Cartesian coordinates, power balance, and shunt
     * admittance to ground
     *
     * A port to the bus reads its voltage and contributes to its power
     * balance. Each balance row is the power into the bus from its shunt and
     * from its ports.
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
      using InternalVariablesT  = BusInternalVariables;
      using ExternalVariablesT  = Utilities::EmptyEnum;
      using InternalEqualitiesT = BusInternalEqualities;
      using ExternalEqualitiesT = Utilities::EmptyEnum;
      using InequalitiesT       = BusInequalities;

      static constexpr std::array<IdxT, 0> NUM_PORT_VARIABLES{};
      static constexpr std::array<IdxT, 0> NUM_PORT_EQUALITIES{};

      explicit Bus(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        initialize() override;

      int evaluateGradient() override;
      int evaluateEqualityJacobian() override;
      int evaluateInequalityJacobian() override;
      int evaluateHessian(RealT sigma, const RealT* lambda_c, const RealT* lambda_d) override;

      [[gnu::always_inline]] void internalEqualities(const ScalarT* x, const ScalarT* external, ScalarT* c) const override;
      [[gnu::always_inline]] void inequalities(const ScalarT* x, const ScalarT* external, ScalarT* d) const override;

      /// Angle reference, with zero imaginary voltage
      bool reference() const
      {
        return reference_;
      }

    private:
      void initializeParameters(const ModelDataT& data);

      RealT Vmin_{ZERO<RealT>};
      RealT Vmax_{UNBOUNDED<RealT>};
      RealT G_{ZERO<RealT>};
      RealT B_{ZERO<RealT>};
      bool  reference_{false};
    };
  } // namespace OptimalPowerFlow
} // namespace GridKit
