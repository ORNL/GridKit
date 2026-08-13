/**
 * @file Genrou.hpp
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Declaration of a GENROU generator model.
 *
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENROU/GenrouData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

// Forward declarations.
namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename scalar_type, typename index_type>
    class BusBase;

    template <typename scalar_type, typename index_type>
    class SignalNode;

    template <typename real_type, typename index_type>
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
      DELTA,  ///< \f$\delta\f$
      OMEGA,  ///< \f$\omega\f$
      EPQ,    ///< \f$E'_q\f$
      PSIPD,  ///< \f$\psi'_d\f$
      PSIPQ,  ///< \f$\psi'_q\f$
      EPD,    ///< \f$E'_d\f$
      PSIPPQ, ///< \f$\psi''_q\f$
      PSIPPD, ///< \f$\psi''_d\f$
      PSIPP,  ///< \f$\psi''\f$
      KSAT,   ///< \f$k_{sat}\f$
      VD,     ///< \f$V_d\f$
      VQ,     ///< \f$V_q\f$
      TE,     ///< \f$T_e\f$
      ID,     ///< \f$I_d\f$
      IQ,     ///< \f$I_q\f$
      IR,     ///< \f$I_r\f$
      II,     ///< \f$I_i\f$
      MAXIMUM,
    };

    /// External variables of a `Genrou`
    enum class GenrouExternalVariables : size_t
    {
      VR,  ///< \f$V_r\f$
      VI,  ///< \f$V_i\f$
      PM,  ///< \f$P_m\f$
      EFD, ///< \f$E_{fd}\f$
      MAXIMUM,
    };

    template <typename scalar_type, typename index_type>
    class Genrou : public Component<scalar_type, index_type>
    {
      using Component<scalar_type, index_type>::gridkit_component_id_;
      using Component<scalar_type, index_type>::alpha_;
      using Component<scalar_type, index_type>::f_;
      using Component<scalar_type, index_type>::nnz_;
      using Component<scalar_type, index_type>::size_;
      using Component<scalar_type, index_type>::tag_;
      using Component<scalar_type, index_type>::abs_tol_;
      using Component<scalar_type, index_type>::time_;
      using Component<scalar_type, index_type>::y_;
      using Component<scalar_type, index_type>::yp_;
      using Component<scalar_type, index_type>::wb_;
      using Component<scalar_type, index_type>::h_;
      using Component<scalar_type, index_type>::J_rows_buffer_;
      using Component<scalar_type, index_type>::J_cols_buffer_;
      using Component<scalar_type, index_type>::J_vals_buffer_;
      using Component<scalar_type, index_type>::freq_system_base_;
      using Component<scalar_type, index_type>::va_system_base_;
      using Component<scalar_type, index_type>::variable_indices_;
      using Component<scalar_type, index_type>::residual_indices_;
      using Component<scalar_type, index_type>::allocated_;

    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename Component<ScalarT, IdxT>::RealT;
      using BusT       = BusBase<ScalarT, IdxT>;
      using ModelDataT = GenrouData<RealT, IdxT>;
      using SignalT    = SignalNode<ScalarT, IdxT>;
      using MonitorT   = Model::VariableMonitor<Genrou, GenrouData>;

      Genrou(BusT* bus);
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
      int setAbsoluteTolerance(RealT) override final;
      int evaluateResidual() override final;

      // Still to be implemented
      int evaluateJacobian() override final;

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
      void initializeParameters(const ModelDataT& data);
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
      __attribute__((always_inline)) inline int evaluateInternalResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);
      __attribute__((always_inline)) inline int evaluateBusResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

    private:
      /* Identification */
      BusT* bus_;
      IdxT  bus_id_{0};

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
