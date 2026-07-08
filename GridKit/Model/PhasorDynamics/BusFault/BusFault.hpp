/* Bus Fault Component - Adam Birchfield */
#pragma once

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/BusFault/BusFaultData.hpp>
#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename scalar_type, typename index_type>
    class BusFault : public Component<scalar_type, index_type>
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
      using Component<scalar_type, index_type>::wb_;
      using Component<scalar_type, index_type>::h_;
      using Component<scalar_type, index_type>::f_;
      using Component<scalar_type, index_type>::J_rows_buffer_;
      using Component<scalar_type, index_type>::J_cols_buffer_;
      using Component<scalar_type, index_type>::J_vals_buffer_;
      using Component<scalar_type, index_type>::variable_indices_;
      using Component<scalar_type, index_type>::residual_indices_;

    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename Component<ScalarT, IdxT>::RealT;
      using BusT       = BusBase<ScalarT, IdxT>;
      using ModelDataT = BusFaultData<RealT, IdxT>;
      using MonitorT   = Model::VariableMonitor<BusFault, BusFaultData>;

      BusFault(BusT* bus);
      BusFault(BusT* bus, RealT R, RealT X, int status);
      BusFault(BusT* bus, const ModelDataT& data);
      ~BusFault();

      int setGridKitComponentID(IdxT) override final;
      int allocate() override final;
      int initialize() override final;
      int tagDifferentiable() override final;
      int setAbsoluteTolerance(RealT rel_tol) override final;
      int evaluateResidual() override final;
      int evaluateJacobian() override final;

      int verify() const override final
      {
        return 0;
      }

      void updateTime(RealT /* t */, RealT /* a */) override final
      {
      }

    public:
      void setR(RealT R)
      {
        R_ = R;
        setDerivedParams();
      }

      void setX(RealT X)
      {
        X_ = X;
        setDerivedParams();
      }

      void setStatus(bool status)
      {
        status_ = status;
      }

      const Model::VariableMonitorBase* getMonitor() const override;

    private:
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
      __attribute__((always_inline)) inline int evaluateBusResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);
      __attribute__((always_inline)) inline int evaluateInternalResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

    private:
      BusT* bus_;
      RealT R_{0.0};
      RealT X_{0.0};
      bool  status_{false};
      IdxT  bus_id_{0};

      /* Derivied parameters */
      RealT B_;
      RealT G_;

      /// Variable monitor
      std::unique_ptr<MonitorT> monitor_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
