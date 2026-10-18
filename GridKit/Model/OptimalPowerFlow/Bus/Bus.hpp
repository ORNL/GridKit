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

    /// Internal constraints of a `Bus`, its power balance first
    enum class BusInternalConstraints : size_t
    {
      P,   ///< Active power balance [p.u.]
      Q,   ///< Reactive power balance [p.u.]
      VM2, ///< \f$V_\mathrm{r}^2 + V_\mathrm{i}^2\f$ Squared voltage magnitude [p.u.]
    };

    /**
     * @brief Bus voltage in Cartesian coordinates, power balance, and shunt
     * admittance to ground
     */
    template <typename scalar_type, typename index_type>
    class Bus final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;
      using ModelDataT = BusData<RealT, IdxT>;

      using ComponentT::g_lower_;
      using ComponentT::g_upper_;
      using ComponentT::tag_;
      using ComponentT::x_;
      using ComponentT::x_lower_;
      using ComponentT::x_upper_;

    public:
      using InternalVariablesT   = BusInternalVariables;
      using ExternalVariablesT   = Utilities::EmptyEnum;
      using InternalConstraintsT = BusInternalConstraints;
      using ExternalConstraintsT = Utilities::EmptyEnum;

      static constexpr std::array<IdxT, 0> NUM_PORT_VARIABLES{};
      static constexpr std::array<IdxT, 0> NUM_PORT_CONSTRAINTS{};

      /// A port to the bus reads its voltage and contributes to its power balance, the constraints before `VM2`
      static constexpr IdxT NUM_SHARED_VARIABLES   = Utilities::enum_size<BusInternalVariables>();
      static constexpr IdxT NUM_SHARED_CONSTRAINTS = static_cast<IdxT>(BusInternalConstraints::VM2);

      explicit Bus(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        tagInequality() override;
      int                        initialize() override;

      int evaluateGradient() override;
      int evaluateJacobian() override;
      int evaluateHessian(RealT sigma) override;

      [[gnu::always_inline]] void internalConstraints(const ScalarT* x, const ScalarT* external, ScalarT* g) const override;

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
