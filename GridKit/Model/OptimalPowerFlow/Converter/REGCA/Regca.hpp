/**
 * @file Regca.hpp
 * @brief Declaration of an optimal power flow REGCA converter.
 */

#pragma once

#include <array>
#include <cstddef>

#include <GridKit/Model/OptimalPowerFlow/Bus/Bus.hpp>
#include <GridKit/Model/OptimalPowerFlow/Component.hpp>
#include <GridKit/Model/OptimalPowerFlow/Converter/REGCA/RegcaData.hpp>
#include <GridKit/Model/OptimalPowerFlow/Signal/Signal.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Converter
    {
      /// Internal variables of a `Regca`, its currents on the component base
      enum class RegcaInternalVariables : size_t
      {
        IP, ///< \f$I_p\f$ Converter active current [p.u.]
        IQ, ///< \f$I_q\f$ Converter reactive current [p.u.]
      };

      /// External variables of a `Regca`, the voltage of its bus, then its signals
      enum class RegcaExternalVariables : size_t
      {
        VR,    ///< \f$V_\mathrm{r}\f$ Real voltage [p.u.]
        VI,    ///< \f$V_\mathrm{i}\f$ Imaginary voltage [p.u.]
        IPCMD, ///< \f$I_p^\mathrm{cmd}\f$ Active-current command [p.u.]
        IQCMD, ///< \f$I_q^\mathrm{cmd}\f$ Reactive-current command [p.u.]
        PBR,   ///< \f$P^\mathrm{br}\f$ Branch active power [p.u.]
        QBR,   ///< \f$Q^\mathrm{br}\f$ Branch reactive power [p.u.]
      };

      /// Internal constraints of a `Regca`
      enum class RegcaInternalConstraints : size_t
      {
        ACTIVE_COMMAND,   ///< Active-current command that the converter current requires [p.u.]
        REACTIVE_COMMAND, ///< Reactive-current command that the converter current requires [p.u.]
        ACTIVE_POWER,     ///< Branch active power [p.u.]
        REACTIVE_POWER,   ///< Branch reactive power [p.u.]
        VT2,              ///< \f$V_T^2\f$ Squared terminal voltage [p.u.]
      };

      /// External constraint contributions of a `Regca`, the power into its bus
      enum class RegcaExternalConstraints : size_t
      {
        DP, ///< \f$\Delta P^\mathrm{bus}\f$ Active power into the bus [p.u.]
        DQ, ///< \f$\Delta Q^\mathrm{bus}\f$ Reactive power into the bus [p.u.]
      };

      /// Renewable converter at steady state: its currents and the commands and power they require
      template <typename scalar_type, typename index_type>
      class Regca final : public Component<scalar_type, index_type>
      {
        using ComponentT = Component<scalar_type, index_type>;
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using RealT      = typename ComponentT::RealT;
        using ModelDataT = RegcaData<RealT, IdxT>;
        using BusT       = Bus<scalar_type, index_type>;
        using SignalT    = Signal<scalar_type, index_type>;

        using ComponentT::g_lower_;
        using ComponentT::g_upper_;
        using ComponentT::tag_;
        using ComponentT::x_upper_;

      public:
        using InternalVariablesT   = RegcaInternalVariables;
        using ExternalVariablesT   = RegcaExternalVariables;
        using InternalConstraintsT = RegcaInternalConstraints;
        using ExternalConstraintsT = RegcaExternalConstraints;

        static constexpr std::array<IdxT, Utilities::enum_size<RegcaPorts>() + Utilities::enum_size<RegcaSignals>()>
            NUM_PORT_VARIABLES{BusT::NUM_SHARED_VARIABLES, SignalT::NUM_SHARED_VARIABLES, SignalT::NUM_SHARED_VARIABLES, SignalT::NUM_SHARED_VARIABLES, SignalT::NUM_SHARED_VARIABLES};
        static constexpr std::array<IdxT, Utilities::enum_size<RegcaPorts>() + Utilities::enum_size<RegcaSignals>()>
            NUM_PORT_CONSTRAINTS{BusT::NUM_SHARED_CONSTRAINTS, 0, 0, 0, 0};

        explicit Regca(const ModelDataT& data);

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

        RealT mva_{ZERO<RealT>};
        bool  sL_{false};
        RealT IL1_{ZERO<RealT>};
        RealT VL1_{ZERO<RealT>};
        RealT VA1_{ZERO<RealT>};
        RealT Vhvmax_{ZERO<RealT>};
      };
    } // namespace Converter
  } // namespace OptimalPowerFlow
} // namespace GridKit
