/**
 * @file Genrou.cpp
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Declaration of a GENROU generator model.
 *
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENROUwS/GenrouData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

// Forward declarations.
namespace GridKit
{
  namespace PhasorDynamics
  {
    template <class ScalarT, typename IdxT>
    class BusBase;

    template <class ScalarT, typename IdxT>
    class SignalNode;

    template <typename RealT, typename IdxT>
    struct GenrouData;
  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Internal variables of a `Genrou`
    enum class GenrouInternalVariables : size_t
    {
      DELTA,  ///< $\delta$
      OMEGA,  ///< $\omega$
      PSIPD,  ///< $\psi'_d$
      PSIPQ,  ///< $\psi'_q$
      EPD,    ///< $E'_d$
      EPQ,    ///< $E'_q$
      VD,     ///< $V_d$
      VQ,     ///< $V_q$
      ID,     ///< $I_d$
      IQ,     ///< $I_q$
      IR,     ///< $I_r$
      II,     ///< $I_i$
      PSIPPQ, ///< $\psi''_q$
      PSIPPD, ///< $\psi''_d$
      PSIPP,  ///< $\psi''$
      TE,     ///< $T_e$
      KSAT,   ///< $k_{sat}$
      MAXIMUM,
    };

    /// External variables of a `Genrou`
    enum class GenrouExternalVariables : size_t
    {
      VR,  ///< $V_r$
      VI,  ///< $V_i$
      PM,  ///< $P_m$
      EFD, ///< $E_{fd}$
      MAXIMUM,
    };

    template <class ScalarT, typename IdxT>
    class Genrou : public Component<ScalarT, IdxT>
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
      using Component<ScalarT, IdxT>::mva_system_base_;
      using Component<ScalarT, IdxT>::variable_indices_;
      using Component<ScalarT, IdxT>::residual_indices_;

    public:
      using RealT           = typename Component<ScalarT, IdxT>::RealT;
      using bus_type        = BusBase<ScalarT, IdxT>;
      using model_data_type = GenrouData<RealT, IdxT>;
      using signal_type     = SignalNode<ScalarT, IdxT>;
      using MonitorT        = Model::VariableMonitor<Genrou, GenrouData>;

      Genrou(bus_type* bus, IdxT unit_id);
      Genrou(bus_type*              bus,
             signal_type*           omega,
             signal_type*           pmech,
             const model_data_type& data);
      Genrou(bus_type*              bus,
             signal_type*           omega,
             signal_type*           pmech,
             signal_type*           efd,
             const model_data_type& data);
      Genrou(bus_type* bus, const model_data_type& data);
      Genrou(bus_type* bus,
             IdxT      unit_id,
             RealT     p0,
             RealT     q0,
             RealT     H,
             RealT     D,
             RealT     Ra,
             RealT     Tdop,
             RealT     Tdopp,
             RealT     Tqopp,
             RealT     Tqop,
             RealT     Xd,
             RealT     Xdp,
             RealT     Xdpp,
             RealT     Xq,
             RealT     Xqp,
             RealT     Xqpp,
             RealT     Xl,
             RealT     S10,
             RealT     S12);
      ~Genrou();

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

      /// Get the `ComponentSignals` from this `Genrou`
      auto getSignals()
          -> ComponentSignals<ScalarT,
                              IdxT,
                              GenrouInternalVariables,
                              GenrouExternalVariables>&
      {
        return signals_;
      }

      const Model::VariableMonitorBase* getMonitor() const override;

    private:
      void initializeParameters(const model_data_type& data);
      /// Associate variable getter functions with enum values
      void initializeMonitor();
      void setDerivedParams();

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
      IdxT      bus_id_{0};
      IdxT      unit_id_; //< @todo this should be removed

      /// Component signal extension
      ComponentSignals<ScalarT, IdxT, GenrouInternalVariables, GenrouExternalVariables> signals_;

      /* Initial terminal conditions */
      RealT p0_{0.0};
      RealT q0_{0.0};

      /* Input parameters */
      RealT H_{0.0};
      RealT D_{0.0};
      RealT Ra_{0.0};
      RealT Tdop_{0.0};
      RealT Tdopp_{0.0};
      RealT Tqopp_{0.0};
      RealT Tqop_{0.0};
      RealT Xd_{0.0};
      RealT Xdp_{0.0};
      RealT Xdpp_{0.0};
      RealT Xq_{0.0};
      RealT Xqp_{0.0};
      RealT Xqpp_{0.0};
      RealT Xl_{0.0};
      RealT S10_{0.0};
      RealT S12_{0.0};
      RealT mva_base_{100.0};

      /* Derivied parameters */
      RealT SA_;
      RealT SB_;
      RealT Xd1_;
      RealT Xd2_;
      RealT Xd3_;
      RealT Xd4_;
      RealT Xd5_;
      RealT Xq1_;
      RealT Xq2_;
      RealT Xq3_;
      RealT Xq4_;
      RealT Xq5_;
      RealT Xqd_;
      RealT G_;
      RealT B_;

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
