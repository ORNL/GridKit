/**
 * @file Machine.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of an EMT synchronous machine.
 */

#pragma once

#include <array>
#include <cstddef>
#include <string_view>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Source/Machine/MachineData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `Machine`, in per unit and of size 1
    enum class MachineInternalVariables : size_t
    {
      theta, ///< \f$\theta\f$ Electrical rotor angle [rad]
      omega, ///< \f$\omega_r\f$ Rotor speed
      psid,  ///< \f$\psi_d\f$ d-axis stator flux linkage
      psiq,  ///< \f$\psi_q\f$ q-axis stator flux linkage
      psi0,  ///< \f$\psi_0\f$ Zero-sequence stator flux linkage
      psifd, ///< \f$\psi_\mathrm{fd}\f$ Field flux linkage
      psi1d, ///< \f$\psi_\mathrm{1d}\f$ d-axis damper flux linkage
      psi1q, ///< \f$\psi_\mathrm{1q}\f$ q-axis damper 1 flux linkage
      psi2q, ///< \f$\psi_\mathrm{2q}\f$ q-axis damper 2 flux linkage
      psiad, ///< \f$\psi_\mathrm{ad}\f$ d-axis air-gap flux linkage
      psiaq, ///< \f$\psi_\mathrm{aq}\f$ q-axis air-gap flux linkage
    };

    /// Outputs of a `Machine`, in per unit and of size 1 unless noted
    enum class MachineOutputs : size_t
    {
      id,      ///< \f$i_d\f$ d-axis stator current
      iq,      ///< \f$i_q\f$ q-axis stator current
      i0,      ///< \f$i_0\f$ Zero-sequence stator current
      ifd,     ///< \f$i_\mathrm{fd}\f$ Field current
      i1d,     ///< \f$i_\mathrm{1d}\f$ d-axis damper current
      i1q,     ///< \f$i_\mathrm{1q}\f$ q-axis damper 1 current
      i2q,     ///< \f$i_\mathrm{2q}\f$ q-axis damper 2 current
      te,      ///< \f$T_e\f$ Electrical torque
      is,      ///< \f$\mathbf{i}_s\f$ Instantaneous stator currents, size 3
      i,       ///< Net currents injected into the bus, phases a, b, c [A]
      p,       ///< Instantaneous three-phase terminal power injected into the grid [W]
      ifd_exc, ///< Field current on the exciter base [p.u.]
    };

    /// Three-phase round-rotor synchronous machine in phase coordinates
    template <typename scalar_type, typename index_type>
    class Machine final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;

      static constexpr size_t PHASES{3};

    public:
      using ModelDataT         = MachineData<RealT, IdxT>;
      using InternalVariablesT = MachineInternalVariables;
      using TerminalsT         = MachineTerminals;
      using InputsT            = MachineInputs;
      using OutputsT           = MachineOutputs;

      static constexpr std::string_view CLASS_NAME = "Machine";

      static constexpr std::array<IdxT, Utilities::enum_size<MachineInternalVariables>()> NUM_INTERNAL_VARIABLES{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
      static constexpr std::array<IdxT, Utilities::enum_size<MachineTerminals>()>         NUM_TERMINAL_VARIABLES{PHASES};
      static constexpr std::array<IdxT, Utilities::enum_size<MachineInputs>()>            NUM_INPUT_VARIABLES{1, 1};
      static constexpr std::array<IdxT, Utilities::enum_size<MachineOutputs>()>           NUM_OUTPUT_VARIABLES{1, 1, 1, 1, 1, 1, 1, 1, PHASES, PHASES, 1, 1};

      explicit Machine(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        tagDifferentiable() override;
      int                        setAbsoluteTolerance(RealT rel_tol) override;

      int evaluateResidual() override;
      int evaluateOutputs() override;
      int evaluateJacobian() override;

      template <typename T>
      [[gnu::always_inline]] void internalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

      template <typename T>
      [[gnu::always_inline]] void externalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

      template <typename T>
      [[gnu::always_inline]] void output(OutputsT output, const T* y, const T* yp, const T* w, const T* wp, T* f) const;

    private:
      // Local indices of the internal variables, in MachineInternalVariables order
      enum : size_t
      {
        THETA,
        OMEGA,
        PSID,
        PSIQ,
        PSI0,
        PSIFD,
        PSI1D,
        PSI1Q,
        PSI2Q,
        PSIAD,
        PSIAQ,
      };

      // Offsets of the external values: terminal `bus`, then inputs `efd` and `pm`
      enum : size_t
      {
        BUS = 0,
        EFD = PHASES,
        PM,
      };

      // Offsets of the winding quantities, in MachineOutputs order
      enum : size_t
      {
        ID,
        IQ,
        I0,
        IFD,
        I1D,
        I1Q,
        I2Q,
        TE,
        IS,
        NUM_WINDING_VALUES = IS + PHASES,
      };

      /// Reciprocal inductances of one rotor axis, zero for invalid data
      struct Axis
      {
        RealT inv_Lm{ZERO<RealT>};
        RealT inv_Ll{ZERO<RealT>};
        RealT inv_L1{ZERO<RealT>};
        RealT inv_L2{ZERO<RealT>};
        RealT inv_leakage_sum{ZERO<RealT>};

        /// Air-gap flux linkage of the axis from its stator and rotor flux linkages
        template <typename T>
        [[gnu::always_inline]] T airGapFlux(T psi, T psi1, T psi2, T inv_ks) const
        {
          return (inv_Ll * psi + inv_L1 * psi1 + inv_L2 * psi2) / (inv_Lm * inv_ks + inv_leakage_sum);
        }
      };

      template <typename T>
      [[gnu::always_inline]] T inverseSaturation(T psi) const;

      template <typename T>
      [[gnu::always_inline]] static std::array<T, PHASES> phaseCurrents(T theta, T id, T iq, T i0);

      template <typename T>
      [[gnu::always_inline]] std::array<T, NUM_WINDING_VALUES> windingValues(const T* y) const;

      RealT S_{ZERO<RealT>};
      RealT V_{ZERO<RealT>};
      RealT f_b_{ZERO<RealT>};
      RealT H_{ZERO<RealT>};
      RealT F_{ZERO<RealT>};
      RealT Rs_{ZERO<RealT>};
      RealT Ll_{ZERO<RealT>};
      RealT Lmd_{ZERO<RealT>};
      RealT Lmq_{ZERO<RealT>};
      RealT L0_{ZERO<RealT>};
      RealT Rfd_{ZERO<RealT>};
      RealT Llfd_{ZERO<RealT>};
      RealT R1d_{ZERO<RealT>};
      RealT Ll1d_{ZERO<RealT>};
      RealT R1q_{ZERO<RealT>};
      RealT Ll1q_{ZERO<RealT>};
      RealT R2q_{ZERO<RealT>};
      RealT Ll2q_{ZERO<RealT>};
      RealT S10_{ZERO<RealT>};
      RealT S12_{ZERO<RealT>};
      RealT Rfe_{ZERO<RealT>};

      RealT omega_b_{ZERO<RealT>};
      RealT Vpk_{ONE<RealT>};
      RealT Ipk_{ZERO<RealT>};
      RealT kfd_{ZERO<RealT>};
      RealT SA_{ZERO<RealT>};
      RealT SB_{ZERO<RealT>};
      RealT Gfe_{ZERO<RealT>}; ///< Iron-loss conductance of a phase [S]

      Axis  d_;
      Axis  q_;
      RealT inv_L0_{ZERO<RealT>};
    };
  } // namespace EMT
} // namespace GridKit
