/**
 * @file Machine.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the EMT synchronous machine model.
 *
 */

#pragma once

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Source/Machine/MachineData.hpp>
#include <GridKit/Model/EMT/ComponentSignals.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

// Forward declarations.
namespace GridKit
{
  namespace EMT
  {
    template <typename real_type, typename index_type>
    struct MachineData;
  } // namespace EMT
} // namespace GridKit

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `Machine`
    enum class MachineInternalVariables : size_t
    {
      THETA, ///< \f$\theta\f$
      OMEGA, ///< \f$\omega_r\f$
      PSID,  ///< \f$\psi_d\f$
      PSIQ,  ///< \f$\psi_q\f$
      PSI0,  ///< \f$\psi_0\f$
      PSIFD, ///< \f$\psi_{fd}\f$
      PSI1D, ///< \f$\psi_{1d}\f$
      PSI1Q, ///< \f$\psi_{1q}\f$
      PSI2Q, ///< \f$\psi_{2q}\f$
      ID,    ///< \f$i_d\f$
      IQ,    ///< \f$i_q\f$
      I0,    ///< \f$i_0\f$
      IFD,   ///< \f$i_{fd}\f$
      I1D,   ///< \f$i_{1d}\f$
      I1Q,   ///< \f$i_{1q}\f$
      I2Q,   ///< \f$i_{2q}\f$
      PSIAD, ///< \f$\psi_{ad}\f$
      PSIAQ, ///< \f$\psi_{aq}\f$
      PSIAT, ///< \f$\psi_{at}\f$
      KS,    ///< \f$K_s\f$
      TE,    ///< \f$T_e\f$
      ISA,   ///< \f$i_{s,a}\f$
      ISB,   ///< \f$i_{s,b}\f$
      ISC,   ///< \f$i_{s,c}\f$
      MAXIMUM,
    };

    /// External variables of a `Machine`
    enum class MachineExternalVariables : size_t
    {
      VA,  ///< \f$v_a\f$
      VB,  ///< \f$v_b\f$
      VC,  ///< \f$v_c\f$
      EFD, ///< \f$E_{fd}\f$
      PM,  ///< \f$P_m\f$
      MAXIMUM,
    };

    /*!
     * @brief Implementation of a three-phase round-rotor synchronous machine.
     *
     * The machine uses fundamental per-unit winding parameters with a field
     * winding, one d-axis damper, and two q-axis dampers. The winding physics
     * rows stay in the rotor dq0 frame; explicit Park residual rows define
     * the dq0 quantities and the instantaneous abc stator currents, and the
     * bus coupling is instantaneous abc SI volts and amps. Round-rotor
     * quadratic saturation applies one factor to both axes.
     */
    template <typename scalar_type, typename index_type>
    class Machine : public Component<scalar_type, index_type>
    {
      using Component<scalar_type, index_type>::gridkit_component_id_;
      using Component<scalar_type, index_type>::size_;
      using Component<scalar_type, index_type>::nnz_;
      using Component<scalar_type, index_type>::time_;
      using Component<scalar_type, index_type>::alpha_;
      using Component<scalar_type, index_type>::y_;
      using Component<scalar_type, index_type>::yp_;
      using Component<scalar_type, index_type>::abs_tol_;
      using Component<scalar_type, index_type>::tag_;
      using Component<scalar_type, index_type>::y_ext_;
      using Component<scalar_type, index_type>::yp_ext_;
      using Component<scalar_type, index_type>::variable_indices_ext_;
      using Component<scalar_type, index_type>::residual_indices_ext_;
      using Component<scalar_type, index_type>::f_ext_;
      using Component<scalar_type, index_type>::f_;
      using Component<scalar_type, index_type>::J_rows_buffer_;
      using Component<scalar_type, index_type>::J_cols_buffer_;
      using Component<scalar_type, index_type>::J_vals_buffer_;
      using Component<scalar_type, index_type>::variable_indices_;
      using Component<scalar_type, index_type>::residual_indices_;
      using Component<scalar_type, index_type>::allocated_;

    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename Component<ScalarT, IdxT>::RealT;
      using ModelDataT = MachineData<RealT, IdxT>;
      using SignalT    = Signal<ScalarT, IdxT>;
      using Port3T     = Port3<ScalarT, IdxT>;
      using MonitorT   = Model::VariableMonitor<Machine, MachineData>;

      Machine();
      Machine(const ModelDataT& data);
      virtual ~Machine();

      virtual int setGridKitComponentID(IdxT) override final;
      virtual int allocate() override final;
      virtual int verify() const override final;
      virtual int initialize() override final;
      virtual int tagDifferentiable() override final;
      virtual int setAbsoluteTolerance(RealT) override final;
      virtual int evaluateInternalResidual() override final;
      virtual int evaluateResidual() override final;
      virtual int evaluateExternalResidual() override final;
      virtual int evaluateJacobian() override final;

      auto getSignals() -> ComponentSignals<ScalarT,
                                            IdxT,
                                            MachineInternalVariables,
                                            MachineExternalVariables>&
      {
        return signals_;
      }

    private:
      void initializeParameters(const ModelDataT& data);
      void initializeMonitor();
      void setDerivedParams();

      ScalarT toMachinePU(ScalarT v) const
      {
        return v / v_peak_base_;
      }

      ScalarT toSystemSI(ScalarT i) const
      {
        return i_peak_base_ * i;
      }

      const Model::VariableMonitorBase* getMonitor() const override;

    protected:
      void gatherExternalVariables() override;

    public:
      __attribute__((always_inline)) inline int evaluateInternalResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);
      __attribute__((always_inline)) inline int evaluateExternalResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

    private:
      /* Initial terminal conditions */
      RealT p0_{0.0};
      RealT q0_{0.0};

      /* Input parameters */
      RealT S_{0.0};
      RealT V_{0.0};
      RealT freq_{0.0};
      RealT H_{0.0};
      RealT Fric_{0.0};
      RealT Rs_{0.0};
      RealT Ll_{0.0};
      RealT Lmd_{0.0};
      RealT Lmq_{0.0};
      RealT L0_{0.0};
      RealT Rfd_{0.0};
      RealT Llfd_{0.0};
      RealT R1d_{0.0};
      RealT Ll1d_{0.0};
      RealT R1q_{0.0};
      RealT Ll1q_{0.0};
      RealT R2q_{0.0};
      RealT Ll2q_{0.0};
      RealT S10_{0.0};
      RealT S12_{0.0};

      /* Derivied parameters */
      RealT omega_base_{0.0};
      RealT v_peak_base_{0.0};
      RealT i_peak_base_{0.0};
      RealT k_fd_{0.0};
      RealT SA_{0.0};
      RealT SB_{0.0};

      /* Setpoints for control variables when no controller is attached */
      ScalarT efd_set_{0.0};
      ScalarT pm_set_{0.0};

      ComponentSignals<ScalarT, IdxT, MachineInternalVariables, MachineExternalVariables> signals_;

      std::unique_ptr<MonitorT> monitor_;
    };

  } // namespace EMT
} // namespace GridKit
