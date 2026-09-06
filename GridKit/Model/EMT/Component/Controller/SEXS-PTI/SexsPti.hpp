/**
 * @file SexsPti.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the SEXS-PTI exciter model.
 */

#pragma once

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/ComponentSignals.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template <typename real_type, typename index_type>
      struct SexsPtiData;
    } // namespace Controller

    template <typename scalar_type, typename index_type>
    class Signal;

  } // namespace EMT
} // namespace GridKit

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      /// Internal variables of a `SexsPti`.
      enum class SexsPtiInternalVariables : size_t
      {
        VR,  ///< Lead-lag block state
        EFD, ///< Exciter field voltage output
        VTR, ///< Terminal voltage error signal
        VTS, ///< Measured terminal voltage magnitude
        MAXIMUM,
      };

      /// External variables of a `SexsPti`.
      enum class SexsPtiExternalVariables : size_t
      {
        VREF, ///< Voltage reference
        VS,   ///< Stabilizer output signal
        VUEL, ///< Under-excitation limiter signal
        VOEL, ///< Over-excitation limiter signal
        VA,   ///< Phase-a terminal voltage in volts
        VB,   ///< Phase-b terminal voltage in volts
        VC,   ///< Phase-c terminal voltage in volts
        MAXIMUM,
      };

      template <typename scalar_type, typename index_type>
      class SexsPti : public Component<scalar_type, index_type>
      {
        using Component<scalar_type, index_type>::gridkit_component_id_;
        using Component<scalar_type, index_type>::alpha_;
        using Component<scalar_type, index_type>::f_;
        using Component<scalar_type, index_type>::nnz_;
        using Component<scalar_type, index_type>::size_;
        using Component<scalar_type, index_type>::tag_;
        using Component<scalar_type, index_type>::abs_tol_;
        using Component<scalar_type, index_type>::time_;
        using Component<scalar_type, index_type>::y_;
        using Component<scalar_type, index_type>::yp_;
        using Component<scalar_type, index_type>::y_ext_;
        using Component<scalar_type, index_type>::yp_ext_;
        using Component<scalar_type, index_type>::variable_indices_ext_;
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
        using ModelDataT = SexsPtiData<RealT, IdxT>;
        using SignalT    = Signal<ScalarT, IdxT>;
        using MonitorT   = Model::VariableMonitor<SexsPti, SexsPtiData>;

        SexsPti();
        explicit SexsPti(const ModelDataT& data);
        ~SexsPti();

        int setGridKitComponentID(IdxT) override final;
        int allocate() override final;
        int verify() const override final;

        int initializationOrder() const noexcept override final
        {
          return 2;
        }

        int initialize() override final;
        int evaluateInternalResidual() override final;
        int tagDifferentiable() override final;
        int setAbsoluteTolerance(RealT rel_tol) override final;
        int evaluateResidual() override final;
        int evaluateJacobian() override final;

        auto getSignals()
            -> ComponentSignals<ScalarT,
                                IdxT,
                                SexsPtiInternalVariables,
                                SexsPtiExternalVariables>&
        {
          return signals_;
        }

        const Model::VariableMonitorBase* getMonitor() const override;

        __attribute__((always_inline)) inline int evaluateInternalResidual(
            const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

      private:
        RealT V_{0};  ///< Rated line-to-line RMS voltage in volts
        RealT Tr_{0}; ///< Optional terminal-voltage measurement lag in seconds

        RealT Ta_{0};
        RealT Tb_{0};
        RealT Te_{0};
        RealT K_{0};
        RealT Efdmax_{0};
        RealT Efdmin_{0};

        int missing_param_count_{0};

        // Runtime connection masks keep the summing junction Enzyme sparse-solvable
        RealT uel_on_{0};
        RealT oel_on_{0};

        ScalarT vref_set_{0};
        ScalarT vs_set_{0};
        ScalarT vuel_set_{0};
        ScalarT voel_set_{0};

        ComponentSignals<ScalarT, IdxT, SexsPtiInternalVariables, SexsPtiExternalVariables> signals_;

        std::unique_ptr<MonitorT> monitor_;

        void initModelParams(const ModelDataT& data);
        void initializeMonitor();

        __attribute__((always_inline)) inline ScalarT voltageMagnitude(const ScalarT* voltage) const
        {
          const ScalarT norm = voltage[0] * voltage[0] + voltage[1] * voltage[1] + voltage[2] * voltage[2];
          return norm == ZERO<RealT> ? ScalarT{ZERO<RealT>} : std::sqrt(norm) / V_;
        }

      protected:
        void gatherExternalVariables() override;
      };

    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
