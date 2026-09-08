/**
 * @file Transformer.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the EMT transformer model.
 *
 */

#pragma once

#include <array>
#include <memory>
#include <string_view>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Transformer/TransformerData.hpp>
#include <GridKit/Model/EMT/ComponentSignals.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

// Forward declarations.
namespace GridKit
{
  namespace EMT
  {
    template <typename real_type, typename index_type>
    struct TransformerData;
  } // namespace EMT
} // namespace GridKit

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `Transformer`
    enum class TransformerInternalVariables : size_t
    {
      I12A,  ///< \f$i_{12,a}\f$
      I12B,  ///< \f$i_{12,b}\f$
      I12C,  ///< \f$i_{12,c}\f$
      PSI1A, ///< \f$\psi_{1,a}\f$
      PSI1B, ///< \f$\psi_{1,b}\f$
      PSI1C, ///< \f$\psi_{1,c}\f$
      PSI2A, ///< \f$\psi_{2,a}\f$
      PSI2B, ///< \f$\psi_{2,b}\f$
      PSI2C, ///< \f$\psi_{2,c}\f$
      E1A,   ///< \f$e_{1,a}\f$
      E1B,   ///< \f$e_{1,b}\f$
      E1C,   ///< \f$e_{1,c}\f$
      E2A,   ///< \f$e_{2,a}\f$
      E2B,   ///< \f$e_{2,b}\f$
      E2C,   ///< \f$e_{2,c}\f$
      IM1A,  ///< \f$i_{m1,a}\f$
      IM1B,  ///< \f$i_{m1,b}\f$
      IM1C,  ///< \f$i_{m1,c}\f$
      IM2A,  ///< \f$i_{m2,a}\f$
      IM2B,  ///< \f$i_{m2,b}\f$
      IM2C,  ///< \f$i_{m2,c}\f$
      IW1A,  ///< \f$i_{w1,a}\f$
      IW1B,  ///< \f$i_{w1,b}\f$
      IW1C,  ///< \f$i_{w1,c}\f$
      IW2A,  ///< \f$i_{w2,a}\f$
      IW2B,  ///< \f$i_{w2,b}\f$
      IW2C,  ///< \f$i_{w2,c}\f$
      MAXIMUM,
    };

    /// External variables of a `Transformer`
    enum class TransformerExternalVariables : size_t
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
     * @brief Implementation of a three-phase bank of two-winding transformers.
     *
     * Each phase is an independent duality-derived pi circuit in transformer
     * per unit: winding resistances outside the ideal ratios, a series
     * leakage reactance between two magnetizing nodes, and a two-slope
     * magnetizing characteristic split between the nodes. Terminal
     * connection maps apply the winding configuration, and the bus coupling
     * is instantaneous abc SI volts and amps.
     */
    template <typename scalar_type, typename index_type>
    class Transformer : public Component<scalar_type, index_type>
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

    public:
      using ScalarT      = scalar_type;
      using IdxT         = index_type;
      using RealT        = typename Component<ScalarT, IdxT>::RealT;
      using ModelDataT   = TransformerData<RealT, IdxT>;
      using Outputs      = typename ModelDataT::Outputs;
      using SignalT      = Signal<ScalarT, IdxT>;
      using MonitorT     = Model::VariableMonitor<Transformer, TransformerData>;
      using PhaseSignals = std::array<SignalT*, 3>;

      Transformer();
      Transformer(const ModelDataT& data);
      virtual ~Transformer();

      void attachTerminal(size_t end, PhaseSignals voltage);

      SignalT& inputSignal(TransformerInputs input)
      {
        return *signals_.getAttachedSignal(static_cast<TransformerExternalVariables>(input));
      }

      SignalT& currentSignal(size_t end, size_t phase)
      {
        return current_.at(end).at(phase);
      }

      virtual int setGridKitComponentID(IdxT) override final;
      virtual int allocate() override final;
      virtual int verify() const override final;

      int initialize(const std::array<RealT, 6>& flux = {});

      typename Component<ScalarT, IdxT>::InitializationPortsT initializationPorts() override
      {
        return {};
      }

      int  initializeState(const std::map<std::string, RealT>& values) override;
      void validateInitialState(const std::map<std::string, RealT>& values) const override;

      /// Initialize from the attached sinusoidal bus-voltage samples.
      int         initializeSteadyState(RealT omega);
      void        assignOutput(Outputs output, SignalT* signal);
      virtual int setAbsoluteTolerance(RealT) override final;
      virtual int evaluateInternalResidual() override final;
      virtual int evaluateResidual() override final;
      virtual int assembleJacobian(RealT y_scale, RealT yp_scale) override final;

      auto getSignals() -> ComponentSignals<ScalarT,
                                            IdxT,
                                            TransformerInternalVariables,
                                            TransformerExternalVariables>&
      {
        return signals_;
      }

    private:
      void initializeParameters(const ModelDataT& data);
      void initializeMonitor();
      void setDerivedParams();

      const ABCMatrix<RealT>& connectionMap(size_t end) const
      {
        if (end == 0)
        {
          return P1_;
        }
        return P2_;
      }

      ScalarT terminalCurrent(size_t end, size_t phase);

      const Model::VariableMonitorBase* getMonitor() const override;

    public:
      __attribute__((always_inline)) inline int evaluateInternalResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

    private:
      /// Two-slope magnetizing characteristic in transformer per unit.
      __attribute__((always_inline)) inline ScalarT magnetizingCurrent(ScalarT psi) const;

      static constexpr std::array<std::string_view, 6> flux_keys_{"psi1a", "psi1b", "psi1c", "psi2a", "psi2b", "psi2c"};

      std::array<std::array<SignalT, 3>, 2> current_;

      /* Input parameters */
      RealT            S_{0.0};
      RealT            V1_{0.0};
      RealT            V2_{0.0};
      RealT            freq_{0.0};
      ABCMatrix<RealT> P1_{{{{1.0, 0.0, 0.0}}, {{0.0, 1.0, 0.0}}, {{0.0, 0.0, 1.0}}}};
      ABCMatrix<RealT> P2_{{{{1.0, 0.0, 0.0}}, {{0.0, 1.0, 0.0}}, {{0.0, 0.0, 1.0}}}};
      RealT            tap_{1.0};
      RealT            R_{0.0};
      RealT            X_{0.0};
      RealT            I0_{0.0};
      RealT            P0_{0.0};
      RealT            knee_{1.2};
      RealT            Lsat_{0.0};
      RealT            split_{0.5};

      /* Derived parameters */
      RealT                omega_base_{0.0};
      std::array<RealT, 2> v_peak_base_{{0.0, 0.0}};
      std::array<RealT, 2> i_peak_base_{{0.0, 0.0}};
      RealT                R1_{0.0};
      RealT                R2_{0.0};
      RealT                Gc_{0.0};
      RealT                Lm_{0.0};
      RealT                inv_Lm_{0.0};
      RealT                k_sat_{0.0};
      std::array<RealT, 2> beta_{{0.0, 0.0}};

      ComponentSignals<ScalarT, IdxT, TransformerInternalVariables, TransformerExternalVariables> signals_;

      std::unique_ptr<MonitorT> monitor_;
    };

  } // namespace EMT
} // namespace GridKit
