/**
 * @file ForcedOscillation.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the forced-oscillation signal source.
 */

#pragma once

#include <cstddef>
#include <memory>

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/PhasorDynamics/SignalSource/ForcedOscillationData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Published signal storage owned by `ForcedOscillation`.
    enum class ForcedOscillationInternalVariables : size_t
    {
      OUTPUT, ///< \f$s_{\mathrm{FO}}\f$ Stateless waveform output
      MAXIMUM ///< Number of published signal values
    };

    /// External signal variables read by `ForcedOscillation`.
    enum class ForcedOscillationExternalVariables : size_t
    {
      MAXIMUM ///< Number of external variables; `ForcedOscillation` has none
    };

    /**
     * @brief Stateless source for a windowed, chirped, and decaying periodic waveform.
     *
     * The source publishes a value owned outside the system DAE. Its signal-node
     * index is always `INVALID_INDEX`, and `updateTime()` refreshes the published
     * value directly from the configured waveform.
     *
     * @tparam scalar_type Plain real or differentiable scalar type.
     * @tparam index_type Integer index type.
     */
    template <typename scalar_type, typename index_type>
    class ForcedOscillation : public Component<scalar_type, index_type>
    {
      using Component<scalar_type, index_type>::allocated_;
      using Component<scalar_type, index_type>::gridkit_component_id_;
      using Component<scalar_type, index_type>::size_;
      using Component<scalar_type, index_type>::time_;

    public:
      using ScalarT            = scalar_type;
      using IdxT               = index_type;
      using RealT              = typename Component<ScalarT, IdxT>::RealT;
      using ModelDataT         = ForcedOscillationData<RealT, IdxT>;
      using MonitorT           = Model::VariableMonitor<ForcedOscillation, ForcedOscillationData>;
      using InternalVariablesT = ForcedOscillationInternalVariables;
      using ExternalVariablesT = ForcedOscillationExternalVariables;

      ForcedOscillation();
      explicit ForcedOscillation(const ModelDataT& data);
      ~ForcedOscillation();

      int  setGridKitComponentID(IdxT component_id) override final;
      int  allocate() override final;
      int  verify() const override final;
      int  initialize() override final;
      int  tagDifferentiable() override final;
      int  setAbsoluteTolerance(RealT rel_tol) override final;
      void updateTime(RealT time, RealT alpha) override final;
      int  evaluateResidual() override final;
      int  evaluateJacobian() override final;

      auto getSignals()
          -> ComponentSignals<ScalarT,
                              IdxT,
                              ForcedOscillationInternalVariables,
                              ForcedOscillationExternalVariables>&;

      const Model::VariableMonitorBase* getMonitor() const override;

    private:
      void  loadRealParameter(const ModelDataT&           data,
                              ForcedOscillationParameters parameter,
                              RealT&                      target,
                              const char*                 name);
      void  initializeParameters(const ModelDataT& data);
      void  initializeMonitor();
      RealT evaluateCarrier(RealT phase) const;
      void  refreshOutput(RealT time);

      RealT                     A_{static_cast<RealT>(0.0)};
      RealT                     f_{static_cast<RealT>(0.0)};
      RealT                     Kf_{static_cast<RealT>(0.0)};
      RealT                     Phi_{static_cast<RealT>(0.0)};
      RealT                     Ton_{static_cast<RealT>(0.0)};
      RealT                     Toff_{static_cast<RealT>(-1.0)};
      RealT                     Tr_{static_cast<RealT>(0.0)};
      RealT                     Tf_{static_cast<RealT>(0.0)};
      RealT                     Kd_{static_cast<RealT>(0.0)};
      ForcedOscillationWaveform waveform_{ForcedOscillationWaveform::SINE};

      IdxT parameter_error_count_{0};

      ScalarT output_{};
      ScalarT envelope_{};
      ScalarT active_{};
      IdxT    output_index_{INVALID_INDEX<IdxT>};

      ComponentSignals<ScalarT,
                       IdxT,
                       ForcedOscillationInternalVariables,
                       ForcedOscillationExternalVariables>
                                signals_;
      std::unique_ptr<MonitorT> monitor_;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
