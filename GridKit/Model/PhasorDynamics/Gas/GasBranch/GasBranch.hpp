/**
 * @file GasBranch.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the gas branch (pipe) model.
 */
#pragma once

#include <memory>

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/Gas/GasBranch/GasBranchData.hpp>
#include <GridKit/Model/PhasorDynamics/Gas/GasBusBase.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Gas
    {
      /**
       * @brief Isothermal gas pipe between two gas buses, divided into N
       * equal segments.
       *
       * The local layout is the interior node pressures p_1 ... p_{N-1}
       * followed by the segment mass fluxes phi_1 ... phi_N. Terminal mass
       * flows are oriented entering the buses.
       */
      template <typename scalar_type, typename index_type>
      class GasBranch : public Component<scalar_type, index_type>
      {
        using Component<scalar_type, index_type>::abs_tol_;
        using Component<scalar_type, index_type>::allocated_;
        using Component<scalar_type, index_type>::alpha_;
        using Component<scalar_type, index_type>::f_;
        using Component<scalar_type, index_type>::gridkit_component_id_;
        using Component<scalar_type, index_type>::h_;
        using Component<scalar_type, index_type>::J_cols_buffer_;
        using Component<scalar_type, index_type>::J_rows_buffer_;
        using Component<scalar_type, index_type>::J_vals_buffer_;
        using Component<scalar_type, index_type>::nnz_;
        using Component<scalar_type, index_type>::residual_indices_;
        using Component<scalar_type, index_type>::size_;
        using Component<scalar_type, index_type>::tag_;
        using Component<scalar_type, index_type>::variable_indices_;
        using Component<scalar_type, index_type>::wb_;
        using Component<scalar_type, index_type>::y_;
        using Component<scalar_type, index_type>::yp_;

      private:
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using RealT      = typename Component<ScalarT, IdxT>::RealT;
        using VectorT    = typename Component<ScalarT, IdxT>::VectorT;
        using GasBusT    = GasBusBase<ScalarT, IdxT>;
        using ModelDataT = GasBranchData<RealT, IdxT>;
        using MonitorT   = Model::VariableMonitor<GasBranch, GasBranchData>;

      public:
        GasBranch(GasBusT* bus1, GasBusT* bus2, const ModelDataT& data);
        ~GasBranch();

        int setGridKitComponentID(IdxT component_id) override final;
        int allocate() override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int setAbsoluteTolerance(RealT rel_tol) override final;
        int evaluateResidual() override final;
        int evaluateJacobian() override final;

        Model::ConfigurationChecks verify() const override final;

        const Model::VariableMonitorBase* getMonitor() const override;

        int evaluateInternalResidual(const ScalarT* y, const ScalarT* yp, const ScalarT* wb, ScalarT* f);
        int evaluateBusResidual(const ScalarT* y, const ScalarT* wb, const ScalarT* wbp, ScalarT* h);

      private:
        void initializeParameters(const ModelDataT& data);
        void initializeMonitor();
        void setDerivedParameters();
        void gatherBusValues();

        ScalarT density(ScalarT p) const;
        ScalarT densitySlope(ScalarT p) const;
        ScalarT meanDensity(ScalarT p_from, ScalarT p_to) const;
        ScalarT potential(ScalarT p) const;
        ScalarT pressureAtPotential(ScalarT target) const;

        IdxT    pIndex(IdxT j) const;
        IdxT    phiIndex(IdxT k) const;
        ScalarT nodePressure(IdxT j, const ScalarT* y, const ScalarT* wb) const;
        IdxT    nodeColumn(IdxT j);
        IdxT    busRow(GasBusT* bus);
        IdxT    busColumn(GasBusT* bus);
        void    addEntry(IdxT row, IdxT column, RealT value);

        GasBusT* bus1_{nullptr};
        GasBusT* bus2_{nullptr};

        RealT L_{ZERO<RealT>};
        RealT D_{ZERO<RealT>};
        RealT lambda_{ZERO<RealT>};
        RealT a_{ZERO<RealT>};
        RealT b1_{ONE<RealT>};
        RealT b2_{ZERO<RealT>};
        RealT dxmax_{static_cast<RealT>(1000.0)};

        IdxT  N_{1};
        RealT dx_{ZERO<RealT>};
        RealT A_{ZERO<RealT>};

        /// Bus pressure derivatives seen by the bus residual
        VectorT wbp_;

        /// Variable monitor
        std::unique_ptr<MonitorT> monitor_;
      };
    } // namespace Gas
  } // namespace PhasorDynamics
} // namespace GridKit
