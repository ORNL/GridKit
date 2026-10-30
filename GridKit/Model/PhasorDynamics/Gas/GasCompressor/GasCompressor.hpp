/**
 * @file GasCompressor.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the gas compressor model.
 */
#pragma once

#include <memory>

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/Gas/GasBusBase.hpp>
#include <GridKit/Model/PhasorDynamics/Gas/GasCompressor/GasCompressorData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Gas
    {
      /**
       * @brief Compressor holding the bus-2 pressure at a fixed ratio to the
       * bus-1 pressure.
       *
       * The compressor has no volume. Its mass flow from bus 1 to bus 2 is
       * algebraic, and the time derivative of the ratio is enforced so that
       * the system stays index 1 when both buses store gas.
       */
      template <typename scalar_type, typename index_type>
      class GasCompressor : public Component<scalar_type, index_type>
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
        using Component<scalar_type, index_type>::y_;
        using Component<scalar_type, index_type>::yp_;

      private:
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using RealT      = typename Component<ScalarT, IdxT>::RealT;
        using VectorT    = typename Component<ScalarT, IdxT>::VectorT;
        using GasBusT    = GasBusBase<ScalarT, IdxT>;
        using ModelDataT = GasCompressorData<RealT, IdxT>;
        using MonitorT   = Model::VariableMonitor<GasCompressor, GasCompressorData>;

      public:
        GasCompressor(GasBusT* bus1, GasBusT* bus2, const ModelDataT& data);
        ~GasCompressor();

        int setGridKitComponentID(IdxT component_id) override final;
        int allocate() override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int setAbsoluteTolerance(RealT rel_tol) override final;
        int evaluateResidual() override final;
        int evaluateJacobian() override final;

        Model::ConfigurationChecks verify() const override final;

        const Model::VariableMonitorBase* getMonitor() const override;

        int evaluateInternalResidual(const ScalarT* wbp, ScalarT* f);
        int evaluateBusResidual(const ScalarT* y, ScalarT* h);

      private:
        void initializeParameters(const ModelDataT& data);
        void gatherBusValues();

        IdxT busRow(GasBusT* bus);
        IdxT busColumn(GasBusT* bus);
        void addEntry(IdxT row, IdxT column, RealT value);

        GasBusT* bus1_{nullptr};
        GasBusT* bus2_{nullptr};

        RealT ratio_{ZERO<RealT>};

        /// Bus pressure derivatives seen by the internal residual
        VectorT wbp_;

        /// Variable monitor
        std::unique_ptr<MonitorT> monitor_;
      };
    } // namespace Gas
  } // namespace PhasorDynamics
} // namespace GridKit
