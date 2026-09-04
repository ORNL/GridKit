/**
 * @file LineLumped.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the EMT lumped line model.
 *
 */

#pragma once

#include <optional>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Line/LineLumped/LineLumpedData.hpp>
#include <GridKit/Model/EMT/ComponentSignals.hpp>
#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFit.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

// Forward declarations.
namespace GridKit
{
  namespace EMT
  {
    template <typename real_type, typename index_type>
    struct LineLumpedData;
  } // namespace EMT
} // namespace GridKit

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `LineLumped`
    enum class LineLumpedInternalVariables : size_t
    {
      I12A,  ///< \f$i_{12,a}\f$
      I12B,  ///< \f$i_{12,b}\f$
      I12C,  ///< \f$i_{12,c}\f$
      ISH1A, ///< \f$i^\mathrm{sh}_{1,a}\f$
      ISH1B, ///< \f$i^\mathrm{sh}_{1,b}\f$
      ISH1C, ///< \f$i^\mathrm{sh}_{1,c}\f$
      ISH2A, ///< \f$i^\mathrm{sh}_{2,a}\f$
      ISH2B, ///< \f$i^\mathrm{sh}_{2,b}\f$
      ISH2C, ///< \f$i^\mathrm{sh}_{2,c}\f$
      MAXIMUM,
    };

    /// External variables of a `LineLumped`
    enum class LineLumpedExternalVariables : size_t
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
     * @brief Implementation of a three-phase lumped EMT line.
     *
     * The initial three-phase formulation uses resistance, inductance,
     * conductance, and capacitance matrices per unit length, scaled by the
     * segment length, with the shunt admittance split between the two
     * terminals. The shunt rows read the terminal voltage derivatives, so
     * connected bus voltages are classified as differential.
     */
    template <typename scalar_type, typename index_type>
    class LineLumped : public Component<scalar_type, index_type>
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
      using Component<scalar_type, index_type>::own_size_;

    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename Component<ScalarT, IdxT>::RealT;
      using ModelDataT = LineLumpedData<RealT, IdxT>;
      using SignalT    = Signal<ScalarT, IdxT>;
      using Port3T     = Port3<ScalarT, IdxT>;
      using VectorFitT = VectorFit<ScalarT, IdxT>;
      using MonitorT   = Model::VariableMonitor<LineLumped, LineLumpedData>;

      LineLumped();
      LineLumped(const ModelDataT& data);
      virtual ~LineLumped();

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
                                            LineLumpedInternalVariables,
                                            LineLumpedExternalVariables>&
      {
        return signals_;
      }

    private:
      void initializeParameters(const ModelDataT& data);
      void initializeMonitor();
      void setDerivedParams();
      bool hasShuntCapacitance() const;

      const Model::VariableMonitorBase* getMonitor() const override;

    public:
      __attribute__((always_inline)) inline int evaluateInternalResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);
      __attribute__((always_inline)) inline int evaluateExternalResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

    private:
      /* Input parameters */
      ABCVector<IdxT>  conductors_{{1, 2, 3}};
      RealT            dx_{0.0};
      ABCMatrix<RealT> Rp_{{{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}}};
      ABCMatrix<RealT> Lp_{{{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}}};
      ABCMatrix<RealT> Gp_{{{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}}};
      ABCMatrix<RealT> Cp_{{{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}}};

      /* Derivied parameters */
      ABCMatrix<RealT> R_{{{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}}};
      ABCMatrix<RealT> L_{{{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}}};
      ABCMatrix<RealT> G_{{{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}}};
      ABCMatrix<RealT> C_{{{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}}};

      /* Setpoints for control variables */
      RealT rl_on_{ONE<RealT>};

      /* Rational submodels */
      std::optional<VectorFitT> z_;
      std::optional<VectorFitT> y1_;
      std::optional<VectorFitT> y2_;
      bool                      fit_ez_singular_{false};
      Port3T                    i12_port_{};
      Port3T                    sh1_rows_port_{};
      Port3T                    sh2_rows_port_{};

      ComponentSignals<ScalarT, IdxT, LineLumpedInternalVariables, LineLumpedExternalVariables> signals_;

      std::unique_ptr<MonitorT> monitor_;
    };

  } // namespace EMT
} // namespace GridKit
