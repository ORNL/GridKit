/**
 * @file Genrou.cpp
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Declaration of a GENROU generator model.
 *
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ConnectedElement.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENROUwS/GenrouData.hpp>

// Forward declarations.
namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename ScalarP, typename IdxP>
    class BusBase;

    template <typename ScalarP, typename IdxP>
    class SignalNode;

    template <typename RealP, typename IdxP>
    struct GenrouData;
  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename, typename>
    class Genrou;

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

    template <typename ScalarP, typename IdxP>
    struct ConnectedElementTraits<Genrou<ScalarP, IdxP>>
    {
      using GenrouT = Genrou<ScalarP, IdxP>;

      using ElementT           = GenrouT;
      using ScalarT            = ScalarP;
      using IdxT               = IdxP;
      using RealT              = typename ScalarTraits<ScalarT>::RealT;
      using ModelDataT         = GenrouData<RealT, IdxT>;
      using InternalVariablesT = GenrouInternalVariables;
      using ExternalVariablesT = GenrouExternalVariables;
      using InterfaceT         = Component<ScalarT, IdxT>;
    };

    template <class ScalarP, typename IdxP>
    class Genrou : public ConnectedElement<Genrou<ScalarP, IdxP>>
    {
      using ConnectedElement<Genrou>::gridkit_component_id_;
      using ConnectedElement<Genrou>::alpha_;
      using ConnectedElement<Genrou>::f_;
      using ConnectedElement<Genrou>::nnz_;
      using ConnectedElement<Genrou>::size_;
      using ConnectedElement<Genrou>::tag_;
      using ConnectedElement<Genrou>::time_;
      using ConnectedElement<Genrou>::y_;
      using ConnectedElement<Genrou>::yp_;
      using ConnectedElement<Genrou>::wb_;
      using ConnectedElement<Genrou>::h_;
      using ConnectedElement<Genrou>::J_;
      using ConnectedElement<Genrou>::J_rows_buffer_;
      using ConnectedElement<Genrou>::J_cols_buffer_;
      using ConnectedElement<Genrou>::J_vals_buffer_;
      using ConnectedElement<Genrou>::mva_system_base_;
      using ConnectedElement<Genrou>::variable_indices_;
      using ConnectedElement<Genrou>::residual_indices_;
      using ConnectedElement<Genrou>::monitor_;
      using ConnectedElement<Genrou>::signals_;

    public:
      using ScalarT    = typename ConnectedElement<Genrou>::ScalarT;
      using IdxT       = typename ConnectedElement<Genrou>::IdxT;
      using RealT      = typename ConnectedElement<Genrou>::RealT;
      using BusT       = BusBase<ScalarT, IdxT>;
      using ModelDataT = GenrouData<RealT, IdxT>;
      using MonitorT   = typename ConnectedElement<Genrou>::MonitorT;
      using SignalT    = SignalNode<ScalarT, IdxT>;

      Genrou(BusT* bus, IdxT unit_id);
      Genrou(BusT*             bus,
             SignalT*          omega,
             SignalT*          pmech,
             const ModelDataT& data);
      Genrou(BusT*             bus,
             SignalT*          omega,
             SignalT*          pmech,
             SignalT*          efd,
             const ModelDataT& data);
      Genrou(BusT* bus, const ModelDataT& data);
      Genrou(BusT* bus,
             IdxT  unit_id,
             RealT p0,
             RealT q0,
             RealT H,
             RealT D,
             RealT Ra,
             RealT Tdop,
             RealT Tdopp,
             RealT Tqopp,
             RealT Tqop,
             RealT Xd,
             RealT Xdp,
             RealT Xdpp,
             RealT Xq,
             RealT Xqp,
             RealT Xqpp,
             RealT Xl,
             RealT S10,
             RealT S12);
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

    private:
      void initializeParameters(const ModelDataT& data);
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
      BusT* bus_;
      IdxT  bus_id_{0};
      IdxT  unit_id_; //< @todo this should be removed

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
    };

  } // namespace PhasorDynamics
} // namespace GridKit
