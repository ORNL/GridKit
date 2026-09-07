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
     * Terminal admittance and legacy series impedance are realized by
     * VectorFit submodels.
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
      using Outputs    = typename ModelDataT::Outputs;
      using SignalT    = Signal<ScalarT, IdxT>;
      using VectorFitT = VectorFit<ScalarT, IdxT>;
      using MonitorT   = Model::VariableMonitor<VoltageSource, VoltageSourceData>;

      VoltageSource();
      VoltageSource(const ModelDataT& data);
      virtual ~VoltageSource();

      SignalT& currentSignal(size_t phase)
      {
        return current_.at(phase);
      }

      virtual int setGridKitComponentID(IdxT) override final;
      virtual int allocate() override final;
      virtual int verify() const override final;

      int initialize(const std::map<Outputs, RealT>& outputs = {});

      int initializeState(const std::map<std::string, RealT>& values) override
      {
        return this->initializeOutputs(*this, values);
      }

      void        assignOutput(Outputs output, SignalT* signal);
      /// Initialize from the attached sinusoidal voltage samples.
      int         initializeSteadyState(RealT omega);
      virtual int setAbsoluteTolerance(RealT) override final;
      virtual int evaluateInternalResidual() override final;
      virtual int evaluateResidual() override final;
      virtual int assembleJacobian(RealT y_scale, RealT yp_scale) override final;

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

    private:
      std::array<SignalT, 3> current_;

      /* Input parameters */
      IdxT             n_phases_{3};
      ABCVector<RealT> E_{{0.0, 0.0, 0.0}};
      ABCVector<RealT> phi_{{0.0, 0.0, 0.0}};
      RealT            omega_{0.0};
      ABCMatrix<RealT> Rs_{{{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}}};
      ABCMatrix<RealT> Ls_{{{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}}};

      /// One for the rational-admittance form of the branch rows, zero otherwise.
      RealT fit_on_{ZERO<RealT>};

      /// Rational source admittance operator
      std::optional<VectorFitT> yfit_;

      /// Legacy series impedance, realized as D + sE without memory states
      std::optional<VectorFitT> zfit_;

      /// The rational admittance linear coefficient must be zero, because
      /// the branch voltage is algebraic
      bool fit_ey_nonzero_{false};

      /// Port over the branch voltage variables read by the rational
      /// admittance
      std::array<SignalT, 3> u_port_{};

      ComponentSignals<ScalarT, IdxT, VoltageSourceInternalVariables, VoltageSourceExternalVariables> signals_;

      std::unique_ptr<MonitorT> monitor_;
      size_t                    jacobian_capacity_{0};
    };

  } // namespace EMT
} // namespace GridKit
