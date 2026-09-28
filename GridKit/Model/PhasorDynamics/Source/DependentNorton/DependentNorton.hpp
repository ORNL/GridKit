/**
 * @file DependentNorton.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the Norton source model.
 */

#pragma once

#include <memory>

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/SignalPorts.hpp>
#include <GridKit/Model/PhasorDynamics/Source/DependentNorton/DependentNortonData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Source
    {
      /// External signal variables of a Norton source.
      enum class DependentNortonExternalVariables : size_t
      {
        INR, ///< \f$I_r^\mathrm{N}\f$ Source current, real component on system base [p.u.]
        INI, ///< \f$I_i^\mathrm{N}\f$ Source current, imaginary component on system base [p.u.]
      };

      /**
       * @brief Controlled current source with a parallel admittance.
       *
       * The model contributes current directly to the connected bus and has
       * no internal variables or equations. Current is positive for injection.
       *
       * @tparam scalar_type Plain real or differentiable scalar type.
       * @tparam index_type Integer index type.
       */
      template <typename scalar_type, typename index_type>
      class DependentNorton : public Component<scalar_type, index_type>
      {
        using Component<scalar_type, index_type>::allocated_;
        using Component<scalar_type, index_type>::gridkit_component_id_;
        using Component<scalar_type, index_type>::h_;
        using Component<scalar_type, index_type>::J_cols_buffer_;
        using Component<scalar_type, index_type>::J_rows_buffer_;
        using Component<scalar_type, index_type>::J_vals_buffer_;
        using Component<scalar_type, index_type>::nnz_;
        using Component<scalar_type, index_type>::size_;
        using Component<scalar_type, index_type>::wb_;
        using Component<scalar_type, index_type>::ws_;
        using Component<scalar_type, index_type>::ws_indices_;
        using Component<scalar_type, index_type>::y_;
        using Component<scalar_type, index_type>::yp_;

      public:
        using ScalarT            = scalar_type;
        using IdxT               = index_type;
        using RealT              = typename Component<ScalarT, IdxT>::RealT;
        using BusT               = BusBase<ScalarT, IdxT>;
        using ModelDataT         = DependentNortonData<RealT, IdxT>;
        using SignalPortsT       = SignalPorts<ScalarT, ModelDataT>;
        using MonitorT           = Model::VariableMonitor<DependentNorton, DependentNortonData>;
        using ExternalVariablesT = DependentNortonExternalVariables;

        DependentNorton(BusT* bus);
        DependentNorton(BusT* bus, const ModelDataT& data);
        ~DependentNorton();

        int setGridKitComponentID(IdxT component_id) override final;
        int allocate() override final;
        int verify() const override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int setAbsoluteTolerance(RealT rel_tol) override final;
        int evaluateResidual() override final;
        int evaluateJacobian() override final;

        SignalPortsT& getPorts()
        {
          return ports_;
        }

        const Model::VariableMonitorBase* getMonitor() const override;

        __attribute__((always_inline)) inline int evaluateBusResidual(
            const ScalarT* y, const ScalarT* yp, const ScalarT* wb, const ScalarT* ws, ScalarT* h);

      private:
        void initializeParameters(const ModelDataT& data);
        void initializeMonitor();

        ScalarT& Vr();
        ScalarT& Vi();
        ScalarT  INr() const;
        ScalarT  INi() const;

        BusT* bus_{nullptr};

        // Input parameters
        RealT G_{0};
        RealT B_{0};

        IdxT parameter_error_count_{2};

        SignalPortsT              ports_;
        std::unique_ptr<MonitorT> monitor_;
      };
    } // namespace Source
  } // namespace PhasorDynamics
} // namespace GridKit
