/**
 * @file GasBusInfinite.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the fixed-pressure gas bus model.
 */
#pragma once

#include <memory>

#include <GridKit/Model/PhasorDynamics/Gas/GasBus/GasBusData.hpp>
#include <GridKit/Model/PhasorDynamics/Gas/GasBusBase.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Gas
    {
      /**
       * @brief Gas bus held at a fixed pressure.
       *
       * The bus has no variables, and it supplies whatever mass flow the
       * connected devices draw.
       */
      template <typename scalar_type, typename index_type>
      class GasBusInfinite : public GasBusBase<scalar_type, index_type>
      {
        using GasBusBase<scalar_type, index_type>::allocated_;
        using GasBusBase<scalar_type, index_type>::bus_id_;
        using GasBusBase<scalar_type, index_type>::nnz_;
        using GasBusBase<scalar_type, index_type>::residual_indices_;
        using GasBusBase<scalar_type, index_type>::size_;
        using GasBusBase<scalar_type, index_type>::variable_indices_;

      private:
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using RealT      = typename GasBusBase<ScalarT, IdxT>::RealT;
        using ModelDataT = GasBusData<RealT, IdxT>;
        using MonitorT   = Model::VariableMonitor<GasBusInfinite, GasBusData>;

      public:
        GasBusInfinite(const ModelDataT& data);
        ~GasBusInfinite();

        int allocate() override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int setAbsoluteTolerance(RealT rel_tol) override final;
        int evaluateResidual() override final;
        int evaluateJacobian() override final;

        Model::ConfigurationChecks verify() const override final;

        const Model::VariableMonitorBase* getMonitor() const override;

        ScalarT& p() override
        {
          return p_;
        }

        ScalarT& q() override
        {
          return q_;
        }

      private:
        RealT   p0_{ZERO<RealT>};
        ScalarT p_{ZERO<RealT>};
        ScalarT q_{ZERO<RealT>};

        /// Variable monitor
        std::unique_ptr<MonitorT> monitor_;
      };
    } // namespace Gas
  } // namespace PhasorDynamics
} // namespace GridKit
