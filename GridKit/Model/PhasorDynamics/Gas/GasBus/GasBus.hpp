/**
 * @file GasBus.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the gas bus model.
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
       * @brief Gas bus with a pressure variable and a mass-flow balance.
       *
       * Mass flow entering the bus is positive. Connected devices add their
       * mass flows after the bus resets the balance to its injection.
       */
      template <typename scalar_type, typename index_type>
      class GasBus : public GasBusBase<scalar_type, index_type>
      {
        using GasBusBase<scalar_type, index_type>::abs_tol_;
        using GasBusBase<scalar_type, index_type>::allocated_;
        using GasBusBase<scalar_type, index_type>::bus_id_;
        using GasBusBase<scalar_type, index_type>::f_;
        using GasBusBase<scalar_type, index_type>::residual_indices_;
        using GasBusBase<scalar_type, index_type>::size_;
        using GasBusBase<scalar_type, index_type>::tag_;
        using GasBusBase<scalar_type, index_type>::variable_indices_;
        using GasBusBase<scalar_type, index_type>::y_;
        using GasBusBase<scalar_type, index_type>::yp_;

      private:
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using RealT      = typename GasBusBase<ScalarT, IdxT>::RealT;
        using ModelDataT = GasBusData<RealT, IdxT>;
        using MonitorT   = Model::VariableMonitor<GasBus, GasBusData>;

      public:
        GasBus(const ModelDataT& data);
        ~GasBus();

        int allocate() override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int setAbsoluteTolerance(RealT rel_tol) override final;
        int evaluateResidual() override final;
        int evaluateJacobian() override final;
        int setWithdrawal(RealT withdrawal) override final;

        Model::ConfigurationChecks verify() const override final;

        const Model::VariableMonitorBase* getMonitor() const override;

        ScalarT& p() override
        {
          return y_.getData()[0];
        }

        ScalarT& q() override
        {
          return f_.getData()[0];
        }

      private:
        RealT p0_{ZERO<RealT>};
        RealT injection_{ZERO<RealT>};

        /// Variable monitor
        std::unique_ptr<MonitorT> monitor_;
      };
    } // namespace Gas
  } // namespace PhasorDynamics
} // namespace GridKit
