#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/PhasorDynamics/LoadZIP/LoadZIPData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

// Forward declarations.
namespace GridKit
{
  namespace PhasorDynamics
  {
    template <class ScalarT, typename IdxT>
    class BusBase;

    template <typename RealT, typename IdxT>
    struct LoadZIPData;
  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    /*!
     * @brief Implementation of a ZIP load.
     *
     */
    template <class ScalarT, typename IdxT>
    class LoadZIP : public Component<ScalarT, IdxT>
    {
      using Component<ScalarT, IdxT>::gridkit_component_id_;
      using Component<ScalarT, IdxT>::size_;
      using Component<ScalarT, IdxT>::nnz_;
      using Component<ScalarT, IdxT>::time_;
      using Component<ScalarT, IdxT>::alpha_;
      using Component<ScalarT, IdxT>::y_;
      using Component<ScalarT, IdxT>::yp_;
      using Component<ScalarT, IdxT>::tag_;
      using Component<ScalarT, IdxT>::wb_;
      using Component<ScalarT, IdxT>::h_;
      using Component<ScalarT, IdxT>::J_rows_buffer_;
      using Component<ScalarT, IdxT>::J_cols_buffer_;
      using Component<ScalarT, IdxT>::J_vals_buffer_;
      using Component<ScalarT, IdxT>::variable_indices_;
      using Component<ScalarT, IdxT>::residual_indices_;

    public:
      using RealT           = typename Component<ScalarT, IdxT>::RealT;
      using bus_type        = BusBase<ScalarT, IdxT>;
      using model_data_type = LoadZIPData<RealT, IdxT>;
      using MonitorT        = Model::VariableMonitor<LoadZIP, LoadZIPData>;

      LoadZIP(bus_type* bus);
      LoadZIP(bus_type* bus, RealT P0, RealT Q0, RealT V0, RealT alphaI,
        RealT alphaP);
      LoadZIP(bus_type* bus, const model_data_type& data);
      virtual ~LoadZIP();

      virtual int setGridKitComponentID(IdxT) override final;
      virtual int allocate() override final;
      virtual int initialize() override final;
      virtual int tagDifferentiable() override final;
      virtual int evaluateResidual() override final;
      virtual int evaluateJacobian() override final;

      virtual int verify() const override final
      {
        return 0;
      }

    public:
      void setP0(RealT P0)
      {
        P0_ = P0;
      }

      void setQ0(RealT Q0)
      {
        Q0_ = Q0;
      }

      void setV0(RealT V0)
      {
        V0_ = V0;
      }
      
      void setalphaI(RealT alphaI)
      {
        alphaI_ = alphaI;
      }
      
      void setalphaP(RealT alphaP)
      {
        alphaP_ = alphaP;
      }

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

      const Model::VariableMonitorBase* getMonitor() const override;

    public:
      __attribute__((always_inline)) inline int evaluateBusResidual(ScalarT*, ScalarT*, ScalarT*, ScalarT*);

    private:
      bus_type* bus_{nullptr};
      RealT     P0_{0};
      RealT     Q0_{0};
      RealT     V0_{1.0};
      RealT     alphaI_{0};
      RealT     alphaP_{0};


      std::unique_ptr<MonitorT> monitor_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
