/**
 * @file BusFault.hpp
 * @author Adam Birchfield
 * @author Luke Lowery
 * @brief Declaration of a bus fault to ground model.
 *
 * The model uses Cartesian coordinates.
 *
 */
#pragma once

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

// Forward declarations.
namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename real_type, typename index_type>
    struct BusFaultData;
  }
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// External variables of a `BusFault`.
    enum class BusFaultExternalVariables : size_t
    {
      STATUS, ///< Status signal masking the fault current contribution
      MAXIMUM,
    };

    /**
     * @brief Implementation of a bus fault to ground.
     *
     * The model is implemented in Cartesian coordinates. The status input
     * signal masks the fault current injected into the bus.
     *
     */
    template <typename scalar_type, typename index_type>
    class BusFault : public Component<scalar_type, index_type>
    {
      using Component<scalar_type, index_type>::gridkit_component_id_;
      using Component<scalar_type, index_type>::size_;
      using Component<scalar_type, index_type>::nnz_;
      using Component<scalar_type, index_type>::y_;
      using Component<scalar_type, index_type>::yp_;
      using Component<scalar_type, index_type>::wb_;
      using Component<scalar_type, index_type>::h_;
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
      using ModelDataT = BusFaultData<RealT, IdxT>;
      using MonitorT   = Model::VariableMonitor<BusFault, BusFaultData>;

      BusFault(BusT* bus);
      BusFault(BusT* bus, const ModelDataT& data);
      virtual ~BusFault();

      virtual int setGridKitComponentID(IdxT) override final;
      virtual int allocate() override final;
      virtual int initialize() override final;
      virtual int tagDifferentiable() override final;
      virtual int setAbsoluteTolerance(RealT rel_tol) override final;
      virtual int evaluateResidual() override final;
      virtual int evaluateJacobian() override final;

      virtual int verify() const override final
      {
        return 0;
      }

      /// Get the `ComponentSignals` from this `BusFault`
      auto getSignals()
          -> ComponentSignals<ScalarT, IdxT, NoVariables, BusFaultExternalVariables>&
      {
        return signals_;
      }

      const Model::VariableMonitorBase* getMonitor() const override;

    private:
      void initializeMonitor();
      void setDerivedParams();
      void faultCurrent(ScalarT& Ir, ScalarT& Ii);

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

    private:
      BusT* bus_;
      RealT R_{0.0};
      RealT X_{0.0};
      IdxT  bus_id_{0};

      /* Derived parameters */
      RealT g_{0.0};
      RealT b_{0.0};

      /// Status signal value masking the fault current contribution
      ScalarT status_{0.0};

      /// Component signal extension
      ComponentSignals<ScalarT, IdxT, NoVariables, BusFaultExternalVariables> signals_;

      /// Variable monitor
      std::unique_ptr<MonitorT> monitor_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
