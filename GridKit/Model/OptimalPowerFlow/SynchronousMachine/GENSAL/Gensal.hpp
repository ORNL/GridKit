/**
 * @file Gensal.hpp
 * @brief Declaration of an optimal power flow salient-pole machine.
 */

#pragma once

#include <array>
#include <cstddef>

#include <GridKit/Model/OptimalPowerFlow/Bus/Bus.hpp>
#include <GridKit/Model/OptimalPowerFlow/Component.hpp>
#include <GridKit/Model/OptimalPowerFlow/Signal/Signal.hpp>
#include <GridKit/Model/OptimalPowerFlow/SynchronousMachine/GENSAL/GensalData.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /// Internal variables of a `Gensal`, its terminal current on the machine base
    enum class GensalInternalVariables : size_t
    {
      IR, ///< \f$I_\mathrm{r}\f$ Real current into the bus [p.u.]
      II, ///< \f$I_\mathrm{i}\f$ Imaginary current into the bus [p.u.]
    };

    /// External variables of a `Gensal`, the voltage of its bus, then its signals
    enum class GensalExternalVariables : size_t
    {
      VR,  ///< \f$V_\mathrm{r}\f$ Real voltage [p.u.]
      VI,  ///< \f$V_\mathrm{i}\f$ Imaginary voltage [p.u.]
      PM,  ///< \f$P_\mathrm{m}\f$ Mechanical power [p.u.]
      EFD, ///< \f$E_\mathrm{fd}\f$ Field voltage [p.u.]
    };

    /// Internal constraints of a `Gensal`
    enum class GensalInternalConstraints : size_t
    {
      MECHANICAL, ///< Mechanical power that the terminal current requires [p.u.]
      FIELD,      ///< Field voltage that the terminal current requires [p.u.]
    };

    /// External constraint contributions of a `Gensal`, the power into its bus
    enum class GensalExternalConstraints : size_t
    {
      DP, ///< \f$\Delta P^\mathrm{bus}\f$ Active power into the bus [p.u.]
      DQ, ///< \f$\Delta Q^\mathrm{bus}\f$ Reactive power into the bus [p.u.]
    };

    /// Salient-pole synchronous machine at steady state: the field voltage that its terminal current requires
    template <typename scalar_type, typename index_type>
    class Gensal final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;
      using ModelDataT = GensalData<RealT, IdxT>;
      using BusT       = Bus<scalar_type, index_type>;
      using SignalT    = Signal<scalar_type, index_type>;

      using ComponentT::tag_;
      using ComponentT::x_lower_;
      using ComponentT::x_upper_;

    public:
      using InternalVariablesT   = GensalInternalVariables;
      using ExternalVariablesT   = GensalExternalVariables;
      using InternalConstraintsT = GensalInternalConstraints;
      using ExternalConstraintsT = GensalExternalConstraints;

      static constexpr std::array<IdxT, Utilities::enum_size<GensalPorts>() + Utilities::enum_size<GensalSignals>()>
          NUM_PORT_VARIABLES{BusT::NUM_SHARED_VARIABLES, SignalT::NUM_SHARED_VARIABLES, SignalT::NUM_SHARED_VARIABLES};
      static constexpr std::array<IdxT, Utilities::enum_size<GensalPorts>() + Utilities::enum_size<GensalSignals>()>
          NUM_PORT_CONSTRAINTS{BusT::NUM_SHARED_CONSTRAINTS, 0, 0};

      explicit Gensal(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        tagInequality() override;
      int                        initialize() override;

      int evaluateGradient() override;
      int evaluateJacobian() override;
      int evaluateHessian(RealT sigma) override;

      [[gnu::always_inline]] void internalConstraints(const ScalarT* x, const ScalarT* external, ScalarT* g) const override;
      [[gnu::always_inline]] void externalConstraints(const ScalarT* x, const ScalarT* external, ScalarT* g) const override;

    private:
      void initializeParameters(const ModelDataT& data);
      void setDerivedParameters();

      RealT Ra_{0.0};
      RealT Xd_{2.1};
      RealT Xdp_{0.2};
      RealT Xq_{0.5};
      RealT S10_{0.0};
      RealT S12_{0.0};
      RealT mva_{100.0};

      RealT SA_{ZERO<RealT>};
      RealT SB_{ZERO<RealT>};
    };
  } // namespace OptimalPowerFlow
} // namespace GridKit
