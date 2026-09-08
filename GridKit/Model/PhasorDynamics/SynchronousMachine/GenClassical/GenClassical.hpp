/**
 * @file GenClassical.hpp
 * @author Abdourahman Barry (abdourahman@vt.edu)
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Declaration of a classical generator model.
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
    /// Internal variables of a `GenClassical`
    enum class GenClassicalInternalVariables : size_t
    {
      DELTA, ///< \f$\delta\f$ rotor angle
      OMEGA, ///< \f$\omega\f$ speed deviation
      TE,    ///< \f$T_e\f$ electrical torque
      IR,    ///< \f$I_r\f$ network real current
      II,    ///< \f$I_i\f$ network imaginary current
      MAXIMUM,
    };

    /// External variables of a `GenClassical`
    enum class GenClassicalExternalVariables : size_t
    {
      VR,  ///< \f$V_r\f$ network real voltage
      VI,  ///< \f$V_i\f$ network imaginary voltage
      PM,  ///< \f$P_m\f$ mechanical power
      EFD, ///< \f$E_{fd}\f$ field voltage
      MAXIMUM,
    };

    template <typename scalar_type, typename index_type>
    class GenClassical : public Component<scalar_type, index_type>
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
      using Component<scalar_type, index_type>::ws_;
      using Component<scalar_type, index_type>::ws_indices_;
      using Component<scalar_type, index_type>::h_;
      using Component<scalar_type, index_type>::J_rows_buffer_;
      using Component<scalar_type, index_type>::J_cols_buffer_;
      using Component<scalar_type, index_type>::J_vals_buffer_;
      using Component<scalar_type, index_type>::freq_system_base_;
      using Component<scalar_type, index_type>::variable_indices_;
      using Component<scalar_type, index_type>::residual_indices_;
      using Component<scalar_type, index_type>::allocated_;

    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename Component<ScalarT, IdxT>::RealT;
      using BusT       = BusBase<ScalarT, IdxT>;
      using ModelDataT = GenClassicalData<RealT, IdxT>;
      using MonitorT   = Model::VariableMonitor<GenClassical, GenClassicalData>;

      GenClassical(BusT* bus, const ModelDataT& data);
      ~GenClassical();

      int setGridKitComponentID(IdxT) override final;
      int allocate() override final;
      int verify() const override final;
      int initialize() override final;
      int tagDifferentiable() override final;
      int setAbsoluteTolerance(RealT rel_tol) override final;
      int evaluateResidual() override final;

      // Still to be implemented
      int evaluateJacobian() override final;

      /// Get the `ComponentSignals` from this `GenClassical`
      auto getSignals()
          -> ComponentSignals<ScalarT,
                              IdxT,
                              GenClassicalInternalVariables,
                              GenClassicalExternalVariables>&
      {
        return signals_;
      }

      const Model::VariableMonitorBase* getMonitor() const override;

    private:
      void initializeParameters(const ModelDataT& data);
      /// Associate variable getter functions with enum values
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

    public:
      __attribute__((always_inline)) inline int evaluateInternalResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);
      __attribute__((always_inline)) inline int evaluateBusResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

    private:
      /* Identification */
      BusT* bus_;

      /// Component signal extension
      ComponentSignals<ScalarT, IdxT, GenClassicalInternalVariables, GenClassicalExternalVariables> signals_;

      /* Initial terminal conditions */
      RealT p0_{0.0};
      RealT q0_{0.0};

      /* Input parameters */
      RealT H_{3.0};
      RealT D_{0.0};
      RealT Ra_{0.0};
      RealT Xdp_{0.2};
      RealT mva_base_{100.0};

      /* Derived parameters */
      RealT G_;
      RealT B_;

      /* Setpoints for control variables (determined at initialization) */
      ScalarT pmech_set_{0.0};
      ScalarT efd_set_{0.0};

      /// Variable monitor
      std::unique_ptr<MonitorT> monitor_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
