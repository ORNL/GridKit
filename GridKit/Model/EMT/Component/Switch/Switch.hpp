/**
 * @file Switch.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the EMT switch model.
 *
 */

#pragma once

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Switch/SwitchData.hpp>
#include <GridKit/Model/EMT/ComponentSignals.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

// Forward declarations.
namespace GridKit
{
  namespace EMT
  {
    template <typename real_type, typename index_type>
    struct SwitchData;
  } // namespace EMT
} // namespace GridKit

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `Switch`
    enum class SwitchInternalVariables : size_t
    {
      I12A, ///< \f$i_{12,a}\f$
      I12B, ///< \f$i_{12,b}\f$
      I12C, ///< \f$i_{12,c}\f$
      MAXIMUM,
    };

    /// External variables of a `Switch`
    enum class SwitchExternalVariables : size_t
    {
      V1A, ///< \f$v_{1,a}\f$
      V1B, ///< \f$v_{1,b}\f$
      V1C, ///< \f$v_{1,c}\f$
      V2A, ///< \f$v_{2,a}\f$
      V2B, ///< \f$v_{2,b}\f$
      V2C, ///< \f$v_{2,c}\f$
      MAXIMUM,
    };

    /*!
     * @brief Implementation of an ideal three-phase EMT switch.
     *
     * The ganged Boolean open command operates all phases and is applied
     * outside the differentiated residual as a constant mask, so the open
     * and closed configurations share one residual row structure. Changing
     * the command changes the Jacobian sparsity pattern; the driver must
     * invalidate the assembled Jacobian structure and reinitialize the
     * integrator after toggling it.
     */
    template <typename scalar_type, typename index_type>
    class Switch : public Component<scalar_type, index_type>
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
      using Component<scalar_type, index_type>::yp_ext_;
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
      using ModelDataT = SwitchData<RealT, IdxT>;
      using SignalT    = SignalNode<ScalarT, IdxT>;
      using Port3T     = Port3<ScalarT, IdxT>;
      using MonitorT   = Model::VariableMonitor<Switch, SwitchData>;

      Switch();
      Switch(const ModelDataT& data);
      virtual ~Switch();

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

      auto getSignals() -> ComponentSignals<ScalarT,
                                            IdxT,
                                            SwitchInternalVariables,
                                            SwitchExternalVariables>&
      {
        return signals_;
      }

      /**
       * @brief Apply the ganged open command.
       *
       * After toggling the command, the driver invalidates the system
       * Jacobian structure and reinitializes the integrator.
       */
      void setOpen(bool open)
      {
        open_ = ZERO<RealT>;
        if (open)
        {
          open_ = ONE<RealT>;
        }
        this->resetJacobianStructure();
      }

      bool isOpen() const
      {
        return open_ != ZERO<RealT>;
      }

    private:
      void initializeParameters(const ModelDataT& data);
      void initializeMonitor();

      const Model::VariableMonitorBase* getMonitor() const override;

    public:
      __attribute__((always_inline)) inline int evaluateInternalResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);
      __attribute__((always_inline)) inline int evaluateExternalResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

    private:
      /* Setpoints for control variables */
      RealT open_{ZERO<RealT>};

      ComponentSignals<ScalarT, IdxT, SwitchInternalVariables, SwitchExternalVariables> signals_;

      std::unique_ptr<MonitorT> monitor_;
    };

  } // namespace EMT
} // namespace GridKit
