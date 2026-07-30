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
    template <typename real_type, typename index_type>
    struct LoadZData;
  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// External variables of a `LoadZ`
    enum class LoadZExternalVariables : size_t
    {
      VR, ///< \f$V_r\f$
      VI, ///< \f$V_i\f$
      MAXIMUM,
    };

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
      using ModelDataT = LoadZData<RealT, IdxT>;
      using MonitorT   = Model::VariableMonitor<LoadZ, LoadZData>;

      LoadZ();
      LoadZ(RealT R, RealT X);
      LoadZ(const ModelDataT& data);
      virtual ~LoadZ();

      virtual int setGridKitComponentID(IdxT) override final;
      virtual int allocate() override final;
      virtual int initialize() override final;
      virtual int tagDifferentiable() override final;
      virtual int setAbsoluteTolerance(RealT) override final;
      virtual int evaluateInternalResidual() override final;
      virtual int evaluateExternalResidual() override final;
      virtual int evaluateJacobian() override final;

      virtual int verify() const override final
      {
        int error_count = 0;
        if (!signals_.template isAttached<LoadZExternalVariables::VR>())
        {
          Log::error() << "LoadZ: VR signal is not attached\n";
          ++error_count;
        }
        if (!signals_.template isAttached<LoadZExternalVariables::VI>())
        {
          Log::error() << "LoadZ: VI signal is not attached\n";
          ++error_count;
        }
        return error_count;
      }

      /// Get the `ComponentSignals` from this component
      auto getSignals()
          -> ComponentSignals<ScalarT, IdxT, NoVariables, LoadZExternalVariables>&
      {
        return signals_;
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
      void gatherExternalVariables();
      void setDerivedParams();

      ScalarT Vr() const
      {
        return signals_.template readExternalVariable<LoadZExternalVariables::VR>();
      }

      ScalarT Vi() const
      {
        return signals_.template readExternalVariable<LoadZExternalVariables::VI>();
      }

      const Model::VariableMonitorBase* getMonitor() const override;

    public:
      __attribute__((always_inline)) inline int evaluateExternalResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);
      __attribute__((always_inline)) inline int evaluateInternalResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

    private:
      RealT R_{0.1};
      RealT X_{0.01};

      /* Derivied parameters */
      RealT b_;
      RealT g_;

      std::unique_ptr<MonitorT> monitor_;

      /// Component signals
      ComponentSignals<ScalarT, IdxT, NoVariables, LoadZExternalVariables> signals_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
