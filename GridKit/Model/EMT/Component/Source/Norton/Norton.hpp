/**
 * @file Norton.hpp
 * @brief Declaration of the EMT Norton model.
 */
#pragma once

#include <array>
#include <map>
#include <string>

#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFit.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    class Norton : public Component<scalar_type, index_type>
    {
      using Component<scalar_type, index_type>::J_cols_buffer_;
      using Component<scalar_type, index_type>::J_rows_buffer_;
      using Component<scalar_type, index_type>::J_vals_buffer_;
      using Component<scalar_type, index_type>::abs_tol_;
      using Component<scalar_type, index_type>::allocated_;
      using Component<scalar_type, index_type>::equation_size_;
      using Component<scalar_type, index_type>::f_;
      using Component<scalar_type, index_type>::gridkit_component_id_;
      using Component<scalar_type, index_type>::nnz_;
      using Component<scalar_type, index_type>::residual_indices_;
      using Component<scalar_type, index_type>::size_;
      using Component<scalar_type, index_type>::tag_;
      using Component<scalar_type, index_type>::variable_indices_;
      using Component<scalar_type, index_type>::y_;
      using Component<scalar_type, index_type>::yp_;

    public:
      using ScalarT      = scalar_type;
      using IdxT         = index_type;
      using RealT        = typename Component<ScalarT, IdxT>::RealT;
      using SignalT      = Signal<ScalarT, IdxT>;
      using PhaseSignals = std::array<SignalT*, 3>;
      using YDataT       = VectorFitData<RealT, IdxT>;

      Norton(const YDataT& Y, PhaseSignals voltage, RealT scale);

      SignalT&       outputSignal(size_t phase);
      const SignalT& outputSignal(size_t phase) const;

      int setGridKitComponentID(IdxT id) override;
      int allocate() override;
      int verify() const override;

      int initialize();
      int initializeState(const std::map<std::string, RealT>& values) override;
      int initializeSteadyState(RealT omega);

      typename Component<ScalarT, IdxT>::InitializationPortsT initializationPorts() override;

      int setAbsoluteTolerance(RealT tolerance) override;
      int evaluateInternalResidual() override;
      int evaluateResidual() override;
      int assembleJacobian(RealT y_scale, RealT yp_scale) override;

    private:
      int initializeShunt();

    public:
      __attribute__((always_inline)) inline int evaluateInternalResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

    private:
      PhaseSignals             voltage_;
      std::array<SignalT, 3>   shunt_;
      VectorFit<ScalarT, IdxT> admittance_;
      size_t                   jacobian_capacity_{0};
    };
  } // namespace EMT
} // namespace GridKit
