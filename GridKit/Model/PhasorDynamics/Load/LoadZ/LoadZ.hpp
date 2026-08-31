#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/PhasorDynamics/Load/LoadZ/LoadZData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

// Forward declarations.
namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename scalar_type, typename index_type>
    class BusBase;

    template <typename real_type, typename index_type>
    struct LoadZData;
  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    /*!
     * @brief Implementation of a constant load.
     *
     */
    template <typename scalar_type, typename index_type>
    class LoadZ : public Component<scalar_type, index_type>
    {
      using Component<scalar_type, index_type>::gridkit_component_id_;
      using Component<scalar_type, index_type>::size_;
      using Component<scalar_type, index_type>::nnz_;
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
      using BusT       = BusBase<ScalarT, IdxT>;
      using ModelDataT = LoadZData<RealT, IdxT>;
      using MonitorT   = Model::VariableMonitor<LoadZ, LoadZData>;

      LoadZ(BusT* bus);
      LoadZ(BusT* bus, RealT R, RealT X);
      LoadZ(BusT* bus, const ModelDataT& data);
      virtual ~LoadZ();

      virtual int setGridKitComponentID(IdxT) override final;
      virtual int allocate() override final;
      virtual int initialize() override final;
      virtual int tagDifferentiable() override final;
      virtual int setAbsoluteTolerance(RealT) override final;
      virtual int evaluateResidual() override final;
      virtual int evaluateJacobian() override final;

      virtual int verify() const override final
      {
        return 0;
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

    private:
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

      const Model::VariableMonitorBase* getMonitor() const override;

    private:
      BusT* bus_{nullptr};
      RealT R_{0.1};
      RealT X_{0.01};

      /* Derivied parameters */
      RealT b_;
      RealT g_;

      std::unique_ptr<MonitorT> monitor_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
