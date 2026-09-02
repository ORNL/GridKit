/**
 * @file LoadZ.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the EMT impedance load model.
 *
 */

#pragma once

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Load/LoadZ/LoadZData.hpp>
#include <GridKit/Model/EMT/ComponentSignals.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

// Forward declarations.
namespace GridKit
{
  namespace EMT
  {
    template <typename real_type, typename index_type>
    struct LoadZData;
  } // namespace EMT
} // namespace GridKit

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `LoadZ`
    enum class LoadZInternalVariables : size_t
    {
      IA, ///< \f$i_a\f$
      IB, ///< \f$i_b\f$
      IC, ///< \f$i_c\f$
      MAXIMUM,
    };

    /// External variables of a `LoadZ`
    enum class LoadZExternalVariables : size_t
    {
      VA, ///< \f$v_a\f$
      VB, ///< \f$v_b\f$
      VC, ///< \f$v_c\f$
      MAXIMUM,
    };

    /*!
     * @brief Implementation of a three-phase EMT impedance load.
     *
     * The initial three-phase formulation uses resistance and inductance
     * matrices, with the injected current as a differential variable.
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
      using ModelDataT = LoadZData<RealT, IdxT>;
      using SignalT    = SignalNode<ScalarT, IdxT>;
      using Port3T     = Port3<ScalarT, IdxT>;
      using MonitorT   = Model::VariableMonitor<LoadZ, LoadZData>;

      LoadZ();
      LoadZ(const ModelDataT& data);
      virtual ~LoadZ();

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
                                            LoadZInternalVariables,
                                            LoadZExternalVariables>&
      {
        return signals_;
      }

    private:
      void initializeParameters(const ModelDataT& data);
      void initializeMonitor();
      bool hasInductance() const;

      const Model::VariableMonitorBase* getMonitor() const override;

    public:
      __attribute__((always_inline)) inline int evaluateInternalResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);
      __attribute__((always_inline)) inline int evaluateExternalResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

    private:
      /* Input parameters */
      ABCMatrix<RealT> R_{{{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}}};
      ABCMatrix<RealT> L_{{{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}}};

      ComponentSignals<ScalarT, IdxT, LoadZInternalVariables, LoadZExternalVariables> signals_;

      std::unique_ptr<MonitorT> monitor_;
    };

  } // namespace EMT
} // namespace GridKit
