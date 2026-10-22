/**
 * @file BranchBreakers.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the BranchBreakers phasor-dynamics branch model.
 */

#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <vector>

#include <GridKit/Model/PhasorDynamics/Branch/BranchBreakers/BranchBreakersData.hpp>
#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/SignalPorts.hpp>
#include <GridKit/Model/VariableMonitor.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename scalar_type, typename index_type>
    class BusBase;

    /// Internal variables of a `BranchBreakers`.
    enum class BranchBreakersInternalVariables : size_t
    {
      Z1, ///< \f$z_1\f$ Differential bus-1 breaker closed-state latch [-]
      Z2, ///< \f$z_2\f$ Differential bus-2 breaker closed-state latch [-]
    };

    /// External signal variables read by a `BranchBreakers`.
    enum class BranchBreakersExternalVariables : size_t
    {
      TRIP1,  ///< \f$s_1\f$ Known bus-1 breaker trip command [-]
      CLOSE1, ///< \f$c_1\f$ Known bus-1 breaker close command [-]
      TRIP2,  ///< \f$s_2\f$ Known bus-2 breaker trip command [-]
      CLOSE2, ///< \f$c_2\f$ Known bus-2 breaker close command [-]
    };

    /**
     * @brief Two-terminal pi branch with a circuit breaker at each terminal.
     *
     * Positive current direction is into the buses. The terminal currents are
     * computed outputs: published on output ports, not system variables.
     *
     * @tparam scalar_type Plain real or differentiable scalar type.
     * @tparam index_type Integer index type.
     */
    template <typename scalar_type, typename index_type>
    class BranchBreakers : public Component<scalar_type, index_type>
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
      using Component<scalar_type, index_type>::output_indices_;
      using Component<scalar_type, index_type>::residual_indices_;
      using Component<scalar_type, index_type>::size_;
      using Component<scalar_type, index_type>::tag_;
      using Component<scalar_type, index_type>::variable_indices_;
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
      using ModelDataT         = BranchBreakersData<RealT, IdxT>;
      using SignalPortsT       = SignalPorts<ScalarT, ModelDataT>;
      using MonitorT           = Model::VariableMonitor<BranchBreakers, BranchBreakersData>;
      using InternalVariablesT = BranchBreakersInternalVariables;
      using ExternalVariablesT = BranchBreakersExternalVariables;

      BranchBreakers(BusT* bus1, BusT* bus2, const ModelDataT& data);
      ~BranchBreakers();

      int setGridKitComponentID(IdxT component_id) override final;
      int allocate() override final;
      int initialize() override final;
      int tagDifferentiable() override final;
      int setAbsoluteTolerance(RealT rel_tol) override final;
      int evaluateResidual() override final;
      int evaluateOutputs() override final;
      int evaluateJacobian() override final;

      Model::ConfigurationChecks verify() const override final;

      SignalPortsT& getPorts()
      {
        return ports_;
      }

      const Model::VariableMonitorBase* getMonitor() const override;

      __attribute__((always_inline)) inline int evaluateInternalResidual(
          const ScalarT* y,
          const ScalarT* yp,
          const ScalarT* wb,
          const ScalarT* ws,
          ScalarT*       f);

      __attribute__((always_inline)) inline int evaluateOutputs(
          const ScalarT* y,
          const ScalarT* yp,
          const ScalarT* wb,
          const ScalarT* ws,
          ScalarT*       o);

      __attribute__((always_inline)) inline int evaluateBusResidual11(
          const ScalarT* y, const ScalarT* yp, const ScalarT* wb, ScalarT* h);
      __attribute__((always_inline)) inline int evaluateBusResidual12(
          const ScalarT* y, const ScalarT* yp, const ScalarT* wb, ScalarT* h);
      __attribute__((always_inline)) inline int evaluateBusResidual21(
          const ScalarT* y, const ScalarT* yp, const ScalarT* wb, ScalarT* h);
      __attribute__((always_inline)) inline int evaluateBusResidual22(
          const ScalarT* y, const ScalarT* yp, const ScalarT* wb, ScalarT* h);

    private:
      using CurrentsT = std::array<ScalarT, Utilities::enum_size<BranchBreakersSignalOutputs>()>;

      void      initializeParameters(const ModelDataT& data);
      void      initializeMonitor();
      void      setDerivedParameters();
      CurrentsT terminalCurrents();

      static __attribute__((always_inline)) inline void evaluateAdmittanceBlock(
          const ScalarT  g,
          const ScalarT  b,
          const ScalarT* wb,
          ScalarT*       h);

      static constexpr RealT TIME_CONSTANT_MINIMUM = static_cast<RealT>(1.0e-3);

      BusT* bus1_{nullptr};
      BusT* bus2_{nullptr};

      RealT R_{ZERO<RealT>};
      RealT X_{ZERO<RealT>};
      RealT G_{ZERO<RealT>};
      RealT B_{ZERO<RealT>};
      RealT Gmag_{ZERO<RealT>};
      RealT Bmag_{ZERO<RealT>};
      RealT tap_{ONE<RealT>};
      RealT phase_{ZERO<RealT>};
      RealT Tbrk_{static_cast<RealT>(0.05)};

      // Derived parameters
      RealT Tlatch_{ZERO<RealT>};
      RealT g11_{ZERO<RealT>};
      RealT b11_{ZERO<RealT>};
      RealT g12_{ZERO<RealT>};
      RealT b12_{ZERO<RealT>};
      RealT g21_{ZERO<RealT>};
      RealT b21_{ZERO<RealT>};
      RealT g22_{ZERO<RealT>};
      RealT b22_{ZERO<RealT>};
      RealT g1o_{ZERO<RealT>}; ///< \f$Y_1^\mathrm{o}\f$, bus-2 breaker open
      RealT b1o_{ZERO<RealT>};
      RealT g2o_{ZERO<RealT>}; ///< \f$Y_2^\mathrm{o}\f$, bus-1 breaker open
      RealT b2o_{ZERO<RealT>};

      /// Global bus-voltage indices in wb_ order
      std::vector<IdxT> wb_indices_;

      SignalPortsT              ports_;
      std::unique_ptr<MonitorT> monitor_;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
