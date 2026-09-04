/**
 * @file Bus.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the EMT bus model.
 *
 */

#pragma once

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Bus/BusData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

// Forward declarations.
namespace GridKit
{
  namespace EMT
  {
    template <typename real_type, typename index_type>
    struct BusData;
  } // namespace EMT
} // namespace GridKit

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `Bus`
    enum class BusInternalVariables : size_t
    {
      VA, ///< \f$v_a\f$
      VB, ///< \f$v_b\f$
      VC, ///< \f$v_c\f$
      MAXIMUM,
    };

    /*!
     * @brief Implementation of a three-phase EMT bus.
     *
     * The bus owns the three phase-voltage variables and the three
     * current-balance residual rows. Connected components accumulate their
     * current injections into the residual rows through the voltage port.
     */
    template <typename scalar_type, typename index_type>
    class Bus : public Component<scalar_type, index_type>
    {
      using Component<scalar_type, index_type>::gridkit_component_id_;
      using Component<scalar_type, index_type>::size_;
      using Component<scalar_type, index_type>::nnz_;
      using Component<scalar_type, index_type>::y_;
      using Component<scalar_type, index_type>::yp_;
      using Component<scalar_type, index_type>::f_;
      using Component<scalar_type, index_type>::tag_;
      using Component<scalar_type, index_type>::abs_tol_;
      using Component<scalar_type, index_type>::y_ext_;
      using Component<scalar_type, index_type>::yp_ext_;
      using Component<scalar_type, index_type>::variable_indices_ext_;
      using Component<scalar_type, index_type>::residual_indices_ext_;
      using Component<scalar_type, index_type>::f_ext_;
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
      using ModelDataT = BusData<RealT, IdxT>;
      using Port3T     = Port3<ScalarT, IdxT>;
      using MonitorT   = Model::VariableMonitor<Bus, BusData>;

      Bus();
      Bus(ScalarT va0, ScalarT vb0, ScalarT vc0);
      Bus(const ModelDataT& data);
      virtual ~Bus();

      virtual int setGridKitComponentID(IdxT) override final;
      virtual int allocate() override final;
      virtual int verify() const override final;
      virtual int initialize() override final;
      virtual int tagDifferentiable() override final;
      virtual int setAbsoluteTolerance(RealT) override final;
      virtual int evaluateInternalResidual() override final;
      virtual int evaluateResidual() override final;
      virtual int evaluateExternalResidual() override final;
      virtual int evaluateJacobian() override final;

      /**
       * @brief The bus connection surface.
       *
       * Each phase signal exposes the phase-voltage variable, its derivative,
       * and the current-balance residual row that connected components
       * accumulate their injections into. The signals are bound in allocate().
       */
      Port3T& voltagePort()
      {
        return v_port_;
      }

    private:
      void initializeMonitor();

      const Model::VariableMonitorBase* getMonitor() const override;

    public:
      __attribute__((always_inline)) inline int evaluateInternalResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);
      __attribute__((always_inline)) inline int evaluateExternalResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

    private:
      /* Initial terminal conditions */
      ScalarT va0_{0.0};
      ScalarT vb0_{0.0};
      ScalarT vc0_{0.0};

      Port3T v_port_{};

      std::unique_ptr<MonitorT> monitor_;
    };

  } // namespace EMT
} // namespace GridKit
