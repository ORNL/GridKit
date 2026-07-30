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
    template <typename real_type, typename index_type>
    struct BusFaultData;
  }
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// External variables of a `BusFault`
    enum class BusFaultExternalVariables : size_t
    {
      VR, ///< \f$V_r\f$
      VI, ///< \f$V_i\f$
      MAXIMUM,
    };

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
      using Component<scalar_type, index_type>::y_ext_;
      using Component<scalar_type, index_type>::variable_indices_ext_;
      using Component<scalar_type, index_type>::residual_indices_ext_;
      using Component<scalar_type, index_type>::f_ext_;
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
      using ModelDataT = BusFaultData<RealT, IdxT>;
      using MonitorT   = Model::VariableMonitor<BusFault, BusFaultData>;

      BusFault();
      BusFault(RealT R, RealT X, int status);
      BusFault(const ModelDataT& data);
      ~BusFault();

      int setGridKitComponentID(IdxT) override final;
      int allocate() override final;
      int initialize() override final;
      int tagDifferentiable() override final;
      int setAbsoluteTolerance(RealT rel_tol) override final;
      int evaluateInternalResidual() override final;
      int evaluateExternalResidual() override final;
      int evaluateJacobian() override final;

      int verify() const override final
      {
        int error_count = 0;
        if (!signals_.template isAttached<BusFaultExternalVariables::VR>())
        {
          Log::error() << "BusFault: VR signal is not attached\n";
          ++error_count;
        }
        if (!signals_.template isAttached<BusFaultExternalVariables::VI>())
        {
          Log::error() << "BusFault: VI signal is not attached\n";
          ++error_count;
        }
        return error_count;
      }

      /// Get the `ComponentSignals` from this component
      auto getSignals()
          -> ComponentSignals<ScalarT, IdxT, NoVariables, BusFaultExternalVariables>&
      {
        return signals_;
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
      void gatherExternalVariables();
      void setDerivedParams();

      ScalarT Vr() const
      {
        return signals_.template readExternalVariable<BusFaultExternalVariables::VR>();
      }

      ScalarT Vi() const
      {
        return signals_.template readExternalVariable<BusFaultExternalVariables::VI>();
      }

    public:
      __attribute__((always_inline)) inline int evaluateExternalResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);
      __attribute__((always_inline)) inline int evaluateInternalResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

    private:
      RealT R_{0.0};
      RealT X_{0.0};
      bool  status_{false};
      IdxT  bus_id_{0};

      /* Derivied parameters */
      RealT B_;
      RealT G_;

      /// Component signals
      ComponentSignals<ScalarT, IdxT, NoVariables, BusFaultExternalVariables> signals_;

      /// Variable monitor
      std::unique_ptr<MonitorT> monitor_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
