/**
 * @file GenClassical.hpp
 * @author Abdourahman Barry (abdourahman@vt.edu)
 * @author Slaven Peles (peless@ornl.gov)
 *
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GenClassical/GenClassicalData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

// Forward declarations.
namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename scalar_type, typename index_type>
    class BusBase;

    template <typename real_type, typename index_type>
    struct GenClassicalData;
  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {

    template <typename scalar_type, typename index_type>
    class GenClassical : public Component<scalar_type, index_type>
    {
      using Component<scalar_type, index_type>::gridkit_component_id_;
      using Component<scalar_type, index_type>::alpha_;
      using Component<scalar_type, index_type>::f_;
      using Component<scalar_type, index_type>::nnz_;
      using Component<scalar_type, index_type>::size_;
      using Component<scalar_type, index_type>::tag_;
      using Component<scalar_type, index_type>::time_;
      using Component<scalar_type, index_type>::y_;
      using Component<scalar_type, index_type>::yp_;
      using Component<scalar_type, index_type>::wb_;
      using Component<scalar_type, index_type>::h_;
      using Component<scalar_type, index_type>::J_;
      using Component<scalar_type, index_type>::J_rows_buffer_;
      using Component<scalar_type, index_type>::J_cols_buffer_;
      using Component<scalar_type, index_type>::J_vals_buffer_;
      using Component<scalar_type, index_type>::freq_system_base_;
      using Component<scalar_type, index_type>::va_system_base_;
      using Component<scalar_type, index_type>::variable_indices_;
      using Component<scalar_type, index_type>::residual_indices_;

    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename Component<ScalarT, IdxT>::RealT;
      using BusT       = BusBase<ScalarT, IdxT>;
      using ModelDataT = GenClassicalData<RealT, IdxT>;
      using MonitorT   = Model::VariableMonitor<GenClassical, GenClassicalData>;

      GenClassical(BusT* bus, int unit_id);
      GenClassical(BusT* bus,
                   int   unit_id,
                   RealT p0,
                   RealT q0,
                   RealT H,
                   RealT D,
                   RealT Ra,
                   RealT Xdp);
      GenClassical(BusT* bus, const ModelDataT& data);
      ~GenClassical();

      int setGridKitComponentID(IdxT) override final;
      int allocate() override final;
      int initialize() override final;
      int tagDifferentiable() override final;
      int evaluateResidual() override final;

      int verify() const override final
      {
        return 0;
      }

      // Still to be implemented
      int evaluateJacobian() override final;

      void setPmech(RealT pmech)
      {
        pmech_set_ = pmech;
      }

      void setEp(RealT ep)
      {
        ep_set_ = ep;
      }

      const Model::VariableMonitorBase* getMonitor() const override;

    private:
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
      FORCE_INLINE int evaluateInternalResidual(ScalarT*, ScalarT*, ScalarT*, ScalarT*);
      FORCE_INLINE int evaluateBusResidual(ScalarT*, ScalarT*, ScalarT*, ScalarT*);

    private:
      /* Identification */
      BusT* bus_;
      IdxT  bus_id_{0};
      int   unit_id_; //< @todo this should be removed

      /* Initial terminal conditions */
      RealT p0_{0.0};
      RealT q0_{0.0};

      /* Input parameters */
      RealT H_{0.0};
      RealT D_{0.0};
      RealT Ra_{0.0};
      RealT Xdp_{0.0};
      RealT mva_base_{100.0};

      /* Derivied parameters */
      RealT G_;
      RealT B_;
      RealT va_machine_base_;

      /* Setpoints for control variables (determined at initialization) */
      ScalarT pmech_set_;
      ScalarT ep_set_;

      /// Variable monitor
      std::unique_ptr<MonitorT> monitor_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
