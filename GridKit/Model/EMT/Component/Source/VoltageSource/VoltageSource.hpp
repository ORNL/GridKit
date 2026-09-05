/**
 * @file VoltageSource.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the EMT voltage source model.
 *
 */

#pragma once

#include <optional>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Source/VoltageSource/VoltageSourceData.hpp>
#include <GridKit/Model/EMT/ComponentSignals.hpp>
#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFit.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

// Forward declarations.
namespace GridKit
{
  namespace EMT
  {
    template <typename real_type, typename index_type>
    struct VoltageSourceData;
  } // namespace EMT
} // namespace GridKit

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `VoltageSource`
    enum class VoltageSourceInternalVariables : size_t
    {
      EA, ///< \f$e_a\f$
      EB, ///< \f$e_b\f$
      EC, ///< \f$e_c\f$
      IA, ///< \f$i_a\f$
      IB, ///< \f$i_b\f$
      IC, ///< \f$i_c\f$
      MAXIMUM,
    };

    /// External variables of a `VoltageSource`
    enum class VoltageSourceExternalVariables : size_t
    {
      VA, ///< \f$v_a\f$
      VB, ///< \f$v_b\f$
      VC, ///< \f$v_c\f$
      MAXIMUM,
    };

    /*!
     * @brief Implementation of a three-phase sinusoidal EMT voltage source.
     *
     * The initial three-phase formulation realizes the terminal admittance as
     * a series resistance and inductance, with the injected current as a
     * differential variable.
     */
    template <typename scalar_type, typename index_type>
    class VoltageSource : public Component<scalar_type, index_type>
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
      using Component<scalar_type, index_type>::equation_size_;

    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename Component<ScalarT, IdxT>::RealT;
      using ModelDataT = VoltageSourceData<RealT, IdxT>;
      using SignalT    = Signal<ScalarT, IdxT>;
      using Port3T     = Port3<ScalarT, IdxT>;
      using VectorFitT = VectorFit<ScalarT, IdxT>;
      using MonitorT   = Model::VariableMonitor<VoltageSource, VoltageSourceData>;

      VoltageSource();
      VoltageSource(const ModelDataT& data);
      virtual ~VoltageSource();

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
                                            VoltageSourceInternalVariables,
                                            VoltageSourceExternalVariables>&
      {
        return signals_;
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
      /* Input parameters */
      ABCVector<RealT> E_{{0.0, 0.0, 0.0}};
      ABCVector<RealT> phi_{{0.0, 0.0, 0.0}};
      RealT            omega_{0.0};
      ABCMatrix<RealT> Rs_{{{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}}};
      ABCMatrix<RealT> Ls_{{{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}}};

      /// Masks selecting the series-matrix or rational-admittance form of
      /// the branch rows; exactly one is one
      RealT rl_on_{ONE<RealT>};
      RealT fit_on_{ZERO<RealT>};

      /// Rational source admittance operator
      std::optional<VectorFitT> yfit_;

      /// The rational admittance linear coefficient must be zero, because
      /// the branch voltage is algebraic
      bool fit_ey_nonzero_{false};

      /// Port over the branch voltage variables read by the rational
      /// admittance
      Port3T u_port_{};

      ComponentSignals<ScalarT, IdxT, VoltageSourceInternalVariables, VoltageSourceExternalVariables> signals_;

      std::unique_ptr<MonitorT> monitor_;
    };

  } // namespace EMT
} // namespace GridKit
