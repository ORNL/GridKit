/* Bus Fault Component - Adam Birchfield */
#pragma once

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

// Forward declaration of BusData structure
namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename RealT, typename IdxT>
    struct BusFaultData;
  }
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <class ScalarT, typename IdxT>
    class BusFault : public Component<ScalarT, IdxT>
    {
      using Component<ScalarT, IdxT>::gridkit_component_id_;
      using Component<ScalarT, IdxT>::alpha_;
      using Component<ScalarT, IdxT>::nnz_;
      using Component<ScalarT, IdxT>::size_;
      using Component<ScalarT, IdxT>::tag_;
      using Component<ScalarT, IdxT>::time_;
      using Component<ScalarT, IdxT>::y_;
      using Component<ScalarT, IdxT>::yp_;
      using Component<ScalarT, IdxT>::wb_;
      using Component<ScalarT, IdxT>::h_;
      using Component<ScalarT, IdxT>::J_rows_buffer_;
      using Component<ScalarT, IdxT>::J_cols_buffer_;
      using Component<ScalarT, IdxT>::J_vals_buffer_;

    public:
      using bus_type = BusBase<ScalarT, IdxT>;
      using RealT    = typename Component<ScalarT, IdxT>::RealT;
      using DataT    = BusFaultData<RealT, IdxT>;
      using MonitorT = Model::VariableMonitor<BusFault, BusFaultData>;

      BusFault(bus_type* bus);
      BusFault(bus_type* bus, RealT R, RealT X, int status);
      BusFault(bus_type* bus, const DataT& data);
      ~BusFault();

      int setGridKitComponentID(IdxT) override final;
      int allocate() override final;
      int initialize() override final;
      int tagDifferentiable() override final;
      int evaluateResidual() override final;
      int evaluateJacobian() override final;

      int verify() const override final
      {
        return 0;
      }

      void updateTime(RealT /* t */, RealT /* a */) override final
      {
      }

      void apply(Action action) override final;

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
      __attribute__((always_inline)) inline int evaluateBusResidual(ScalarT*, ScalarT*, ScalarT*, ScalarT*);

    private:
      bus_type* bus_;
      RealT     R_{0.0};
      RealT     X_{0.0};
      bool      status_{false};
      IdxT      bus_id_{0};

      /* Derivied parameters */
      RealT B_;
      RealT G_;

      /// Variable monitor
      std::unique_ptr<MonitorT> monitor_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
