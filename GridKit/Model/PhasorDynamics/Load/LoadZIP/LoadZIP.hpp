#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/PhasorDynamics/Load/LoadZIP/LoadZIPData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

// Forward declarations.
namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename scalar_type, typename index_type>
    class BusBase;

    template <typename real_type, typename index_type>
    struct LoadZIPData;
  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// External variables of a ZIP load
    enum class LoadZIPExternalVariables : size_t
    {
      P,      ///< Initial terminal active-power injection
      Q,      ///< Initial terminal reactive-power injection
      ONLINE, ///< Connection status (zero is disconnected, nonzero is connected)
      MAXIMUM
    };

    /*!
     * @brief Implementation of a ZIP load.
     *
     */
    template <typename scalar_type, typename index_type>
    class LoadZIP : public Component<scalar_type, index_type>
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
      using Component<scalar_type, index_type>::allocated_;

    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename Component<ScalarT, IdxT>::RealT;
      using BusT       = BusBase<ScalarT, IdxT>;
      using ModelDataT = LoadZIPData<RealT, IdxT>;
      using MonitorT   = Model::VariableMonitor<LoadZIP, LoadZIPData>;
      using SignalsT   = ComponentSignals<ScalarT,
                                          IdxT,
                                          NoVariables,
                                          LoadZIPExternalVariables>;

      LoadZIP(BusT* bus);
      LoadZIP(BusT* bus, RealT alphaI, RealT alphaP);
      LoadZIP(BusT* bus, const ModelDataT& data);
      ~LoadZIP();

      int setGridKitComponentID(IdxT) override final;
      int allocate() override final;
      int initialize() override final;
      int tagDifferentiable() override final;
      int setAbsoluteTolerance(RealT rel_tol) override final;
      int evaluateResidual() override final;
      int evaluateJacobian() override final;

      /// Get the signal connections for this ZIP load.
      SignalsT& getSignals()
      {
        return signals_;
      }

      int verify() const override final
      {
        return 0;
      }

    public:
      void setAlphaI(RealT alphaI)
      {
        alphaI_ = alphaI;
        setDerivedParams();
      }

      void setAlphaP(RealT alphaP)
      {
        alphaP_ = alphaP;
        setDerivedParams();
      }

    private:
      void    initializeParameters(const ModelDataT& data);
      void    initializeMonitor();
      void    setDerivedParams();
      void    setInputDispatchAtVoltage(ScalarT vr, ScalarT vi);
      ScalarT online() const;

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
      __attribute__((always_inline)) inline int evaluateBusResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);
      __attribute__((always_inline)) inline int evaluateInternalResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

    private:
      BusT*   bus_{nullptr};
      ScalarT p_{0.0};
      ScalarT q_{0.0};
      /// ZIP anchor voltage, derived from the bus voltage at initialization
      RealT   Vnom_{1.0};
      RealT   alphaI_{0};
      RealT   alphaP_{0};
      ScalarT G_{0.0};
      ScalarT B_{0.0};
      RealT   alphaZ_{1.0};

      SignalsT signals_;

      std::unique_ptr<MonitorT> monitor_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
