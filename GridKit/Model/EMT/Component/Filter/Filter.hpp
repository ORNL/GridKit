/**
 * @file Filter.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Three-phase LCL filter
 *
 */

#pragma once

#include <array>
#include <complex>
#include <memory>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Filter/FilterData.hpp>
#include <GridKit/Model/EMT/ComponentSignals.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `Filter`
    enum class FilterInternalVariables : size_t
    {
      IA,  ///< \f$i_a\f$
      IB,  ///< \f$i_b\f$
      IC,  ///< \f$i_c\f$
      VOA, ///< \f$v_{\mathrm{o},a}\f$
      VOB, ///< \f$v_{\mathrm{o},b}\f$
      VOC, ///< \f$v_{\mathrm{o},c}\f$
      IGA, ///< \f$i_{g,a}\f$
      IGB, ///< \f$i_{g,b}\f$
      IGC, ///< \f$i_{g,c}\f$
      MAXIMUM,
    };

    /// External variables of a `Filter`
    enum class FilterExternalVariables : size_t
    {
      VA, ///< \f$v_a\f$
      VB, ///< \f$v_b\f$
      VC, ///< \f$v_c\f$
      EA, ///< \f$e_a\f$
      EB, ///< \f$e_b\f$
      EC, ///< \f$e_c\f$
      MAXIMUM,
    };

    /**
     * @brief LCL filter with coupled three-phase resistance, inductance and capacitance.
     */
    template <typename scalar_type, typename index_type>
    class Filter : public Component<scalar_type, index_type>
    {
      using Component<scalar_type, index_type>::gridkit_component_id_;
      using Component<scalar_type, index_type>::size_;
      using Component<scalar_type, index_type>::equation_size_;
      using Component<scalar_type, index_type>::nnz_;
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
      using ModelDataT   = FilterData<RealT, IdxT>;
      using Outputs      = typename ModelDataT::Outputs;
      using SignalT      = Signal<ScalarT, IdxT>;
      using MonitorT     = Model::VariableMonitor<Filter, FilterData>;
      using PhaseSignals = std::array<SignalT*, 3>;

      explicit Filter(const ModelDataT& data);
      virtual ~Filter();

      void attachInput(PhaseSignals voltage, PhaseSignals source);
      void assignOutput(Outputs output, SignalT* signal);

      SignalT& outputSignal(Outputs output);

      SignalT& inputSignal(FilterInputs input)
      {
        return *signals_.getAttachedSignal(static_cast<FilterExternalVariables>(input));
      }

      /// Grid-side current, positive into the connected bus.
      SignalT& currentSignal(size_t phase)
      {
        return output_.at(static_cast<size_t>(Outputs::iga) + phase);
      }

      virtual int setGridKitComponentID(IdxT id) override final;
      virtual int allocate() override final;
      virtual int verify() const override final;

      int initialize(const std::map<Outputs, RealT>& outputs = {});

      int  initializeState(const std::map<std::string, RealT>& values) override;
      int  initializeState(const std::map<std::string, RealT>& values, RealT omega) override;
      void validateInitialState(const std::map<std::string, RealT>& values) const override;

      typename Component<ScalarT, IdxT>::InitializationPortsT initializationPorts() override;
      void                                                    prepareInitialization(typename Component<ScalarT, IdxT>::InitialStateT& initial) override;

      virtual int setAbsoluteTolerance(RealT tolerance) override final;
      virtual int evaluateInternalResidual() override final;
      virtual int evaluateResidual() override final;
      virtual int assembleJacobian(RealT y_scale, RealT yp_scale) override final;

      auto getSignals() -> ComponentSignals<ScalarT,
                                            IdxT,
                                            FilterInternalVariables,
                                            FilterExternalVariables>&
      {
        return signals_;
      }

    private:
      std::array<ABCVector<std::complex<RealT>>, 4> operatingPoint(const ABCVector<RealT>& voltage, const ABCVector<RealT>& current, RealT omega) const;
      void                                          initializeParameters(const ModelDataT& data);
      void                                          initializeMonitor();

      const Model::VariableMonitorBase* getMonitor() const override;

    public:
      __attribute__((always_inline)) inline int evaluateInternalResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

    private:
      /* Input parameters */
      ABCMatrix<RealT> Rs_{};
      ABCMatrix<RealT> Ls_{};
      ABCMatrix<RealT> C_{};
      ABCMatrix<RealT> Rg_{};
      ABCMatrix<RealT> Lg_{};

      ComponentSignals<ScalarT, IdxT, FilterInternalVariables, FilterExternalVariables> signals_;

      std::array<SignalT, static_cast<size_t>(FilterInternalVariables::MAXIMUM)>  output_;
      std::array<SignalT*, static_cast<size_t>(FilterInternalVariables::MAXIMUM)> alias_{};
      std::unique_ptr<MonitorT>                                                   monitor_;
    };
  } // namespace EMT
} // namespace GridKit
