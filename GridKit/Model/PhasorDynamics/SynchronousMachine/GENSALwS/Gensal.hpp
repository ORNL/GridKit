/**
 * @file Gensal.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of a GENSAL generator model.
 *
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENSALwS/GensalData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

// Forward declarations.
namespace GridKit
{
  namespace PhasorDynamics
  {
    template <class ScalarT, typename IdxT>
    class BusBase;

    template <typename RealT, typename IdxT>
    struct GensalData;
  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Internal variables of a `Gensal`
    enum class GensalInternalVariables : size_t
    {
      DELTA,  ///< rotor angle
      OMEGA,  ///< speed deviation
      EPQ,    ///< q-axis transient voltage
      PSIPD,  ///< d-axis transient flux
      PSIPPQ, ///< q-axis subtransient flux
      PSIPPD, ///< d-axis subtransient flux
      KSAT,   ///< saturation signal
      VD,     ///< d-axis terminal voltage
      VQ,     ///< q-axis terminal voltage
      TE,     ///< electrical torque
      ID,     ///< d-axis current
      IQ,     ///< q-axis current
      IR,     ///< network real current
      II,     ///< network imaginary current
      INR,    ///< Norton source real current
      INI,    ///< Norton source imaginary current
      MAXIMUM,
    };

    /// External variables of a `Gensal`
    enum class GensalExternalVariables : size_t
    {
      VR,  ///< network real voltage
      VI,  ///< network imaginary voltage
      PM,  ///< mechanical power
      EFD, ///< field voltage
      MAXIMUM,
    };

    template <class ScalarT, typename IdxT>
    class Gensal : public Component<ScalarT, IdxT>
    {
      using Component<ScalarT, IdxT>::gridkit_component_id_;
      using Component<ScalarT, IdxT>::alpha_;
      using Component<ScalarT, IdxT>::f_;
      using Component<ScalarT, IdxT>::nnz_;
      using Component<ScalarT, IdxT>::size_;
      using Component<ScalarT, IdxT>::tag_;
      using Component<ScalarT, IdxT>::time_;
      using Component<ScalarT, IdxT>::y_;
      using Component<ScalarT, IdxT>::yp_;
      using Component<ScalarT, IdxT>::wb_;
      using Component<ScalarT, IdxT>::h_;
      using Component<ScalarT, IdxT>::J_;
      using Component<ScalarT, IdxT>::J_rows_buffer_;
      using Component<ScalarT, IdxT>::J_cols_buffer_;
      using Component<ScalarT, IdxT>::J_vals_buffer_;
      using Component<ScalarT, IdxT>::freq_system_base_;
      using Component<ScalarT, IdxT>::va_system_base_;
      using Component<ScalarT, IdxT>::variable_indices_;
      using Component<ScalarT, IdxT>::residual_indices_;

    public:
      using RealT           = typename Component<ScalarT, IdxT>::RealT;
      using bus_type        = BusBase<ScalarT, IdxT>;
      using model_data_type = GensalData<RealT, IdxT>;
      using MonitorT        = Model::VariableMonitor<Gensal, GensalData>;

      Gensal(bus_type* bus, const model_data_type& data);
      ~Gensal();

      int setGridKitComponentID(IdxT) override final;
      int allocate() override final;
      int verify() const override final;
      int initialize() override final;
      int tagDifferentiable() override final;
      int evaluateResidual() override final;

      // Still to be implemented
      int evaluateJacobian() override final;

      // Temporary access functions for governor
      // Should be abstracted
      ScalarT getSpeed();
      ScalarT getTorque();

      /// Get the `ComponentSignals` from this `Gensal`
      auto getSignals()
          -> ComponentSignals<ScalarT,
                              IdxT,
                              GensalInternalVariables,
                              GensalExternalVariables>&
      {
        return signals_;
      }

      const Model::VariableMonitorBase* getMonitor() const override;

    private:
      void initializeParameters(const model_data_type& data);
      /// Associate variable getter functions with enum values
      void initializeMonitor();
      void setDerivedParams();

      /**
       * @brief Convert per-unit current or power from system base to machine base.
       *
       * @note For terminal-current quantities, this scaling assumes the machine
       * voltage base matches the interfacing bus voltage base. A voltage-base
       * mismatch is not a concern here because the model is formulated at the
       * machine terminals using the connected bus voltage base.
       */
      ScalarT toMachineBase(ScalarT value) const;

      /**
       * @brief Convert per-unit current or power from machine base to system base.
       *
       * @note For terminal-current quantities, this scaling assumes the machine
       * voltage base matches the interfacing bus voltage base. A voltage-base
       * mismatch is not a concern here because the model is formulated at the
       * machine terminals using the connected bus voltage base.
       */
      ScalarT toSystemBase(ScalarT value) const;

      ScalarT& Vr()
      {
        return bus_->Vr();
      }

      ScalarT& Vi()
      {
        return bus_->Vi();
      }

      ScalarT& Ir()
      {
        return bus_->Ir();
      }

      ScalarT& Ii()
      {
        return bus_->Ii();
      }

    public:
      __attribute__((always_inline)) inline int evaluateInternalResidual(ScalarT*, ScalarT*, ScalarT*, ScalarT*, ScalarT*);
      __attribute__((always_inline)) inline int evaluateBusResidual(ScalarT*, ScalarT*, ScalarT*, ScalarT*);

    private:
      /* Identification */
      bus_type* bus_;

      /// Component signal extension
      ComponentSignals<ScalarT, IdxT, GensalInternalVariables, GensalExternalVariables> signals_;

      /* Initial terminal conditions */
      RealT p0_{0.0};
      RealT q0_{0.0};

      /* Input parameters */
      RealT H_{3.0};
      RealT D_{0.0};
      RealT Ra_{0.0};
      RealT Tdop_{7.0};
      RealT Tdopp_{0.04};
      RealT Tqopp_{0.05};
      RealT Xd_{2.1};
      RealT Xdp_{0.2};
      RealT Xdpp_{0.18};
      RealT Xq_{0.5};
      RealT Xl_{0.15};
      RealT S10_{0.0};
      RealT S12_{0.0};
      RealT mva_base_{100.0};

      /* Derived parameters */
      RealT SA_;
      RealT SB_;
      RealT Xd1_;
      RealT Xd2_;
      RealT Xd3_;
      RealT Xd4_;
      RealT Xd5_;
      RealT Xq2_;
      RealT G_;
      RealT B_;
      RealT va_machine_base_;

      /* Setpoints for control variables (determined at initialization) */
      ScalarT pmech_set_{0.0}; // TODO remove default initialization and ensure this gets set
      ScalarT efd_set_{0.0};   // TODO remove default initialization and ensure this gets set

      /* Local copies of signal variables */
      std::vector<ScalarT> ws_;
      std::vector<IdxT>    ws_indices_;

      /// Variable monitor
      std::unique_ptr<MonitorT> monitor_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
