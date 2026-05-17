/**
 * @file Regca.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the REGCA phasor-dynamics converter model.
 */

#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REGCA/RegcaData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <class ScalarT, typename IdxT>
    class BusBase;

    template <class ScalarT, typename IdxT>
    class SignalNode;
  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Converter
    {
      /// Internal variables of a `Regca`
      enum class RegcaInternalVariables : size_t
      {
        VM,      ///< Filtered terminal voltage
        IQ,      ///< Reactive-current state
        IP,      ///< Active-current state
        VT,      ///< Terminal voltage magnitude
        II,      ///< Imaginary injected current
        IQEXTRA, ///< HVRCM extra reactive current
        IL,      ///< LVPL upper-limit current curve
        IR,      ///< Real injected current
        LP,      ///< Active-current lower rate bound
        UP,      ///< Active-current upper rate bound
        MAXIMUM,
      };

      /// External variables of a `Regca`
      enum class RegcaExternalVariables : size_t
      {
        IPCMD, ///< Active-current command signal
        IQCMD, ///< Reactive-current command signal
        MAXIMUM,
      };

      template <class ScalarT, typename IdxT>
      class Regca : public Component<ScalarT, IdxT>
      {
        using Component<ScalarT, IdxT>::gridkit_component_id_;
        using Component<ScalarT, IdxT>::alpha_;
        using Component<ScalarT, IdxT>::f_;
        using Component<ScalarT, IdxT>::h_;
        using Component<ScalarT, IdxT>::J_;
        using Component<ScalarT, IdxT>::J_cols_buffer_;
        using Component<ScalarT, IdxT>::J_rows_buffer_;
        using Component<ScalarT, IdxT>::J_vals_buffer_;
        using Component<ScalarT, IdxT>::mva_system_base_;
        using Component<ScalarT, IdxT>::nnz_;
        using Component<ScalarT, IdxT>::residual_indices_;
        using Component<ScalarT, IdxT>::size_;
        using Component<ScalarT, IdxT>::tag_;
        using Component<ScalarT, IdxT>::time_;
        using Component<ScalarT, IdxT>::variable_indices_;
        using Component<ScalarT, IdxT>::wb_;
        using Component<ScalarT, IdxT>::y_;
        using Component<ScalarT, IdxT>::yp_;

      public:
        using RealT           = typename Component<ScalarT, IdxT>::RealT;
        using bus_type        = BusBase<ScalarT, IdxT>;
        using signal_type     = SignalNode<ScalarT, IdxT>;
        using model_data_type = RegcaData<RealT, IdxT>;
        using MonitorT        = Model::VariableMonitor<Regca, RegcaData>;

        Regca(bus_type* bus);
        Regca(bus_type* bus, const model_data_type& data);
        ~Regca();

        int setGridKitComponentID(IdxT) override final;
        int allocate() override final;
        int verify() const override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int evaluateResidual() override final;
        int evaluateJacobian() override final;

        auto getSignals()
            -> ComponentSignals<ScalarT,
                                IdxT,
                                RegcaInternalVariables,
                                RegcaExternalVariables>&
        {
          return signals_;
        }

        const Model::VariableMonitorBase* getMonitor() const override;

        __attribute__((always_inline)) inline int evaluateInternalResidual(
            ScalarT*, ScalarT*, ScalarT*, ScalarT*, ScalarT*);
        __attribute__((always_inline)) inline int evaluateBusResidual(
            ScalarT*, ScalarT*, ScalarT*, ScalarT*);

      private:
        void initializeParameters(const model_data_type& data);
        void initializeMonitor();

        ScalarT& Vr()
        {
          return bus_->Vr();
        }

        ScalarT& Vi()
        {
          return bus_->Vi();
        }

        ScalarT& Ir()
        {
          return bus_->Ir();
        }

        ScalarT& Ii()
        {
          return bus_->Ii();
        }

        bus_type* bus_{nullptr};

        RealT P0_{0};
        RealT Q0_{0};
        RealT Sconv_{0};
        RealT Tg_{0};
        RealT TM_{0};
        RealT Rqmax_{0};
        RealT Rqmin_{0};
        RealT Rpmax_{0};
        bool  sL_{false};
        RealT IL1_{0};
        RealT VL0_{0};
        RealT VL1_{0};
        RealT VA0_{0};
        RealT VA1_{0};
        RealT Vhvmax_{0};
        IdxT  bus_id_{0};

        ComponentSignals<ScalarT, IdxT, RegcaInternalVariables, RegcaExternalVariables> signals_;
        std::unique_ptr<MonitorT>                                                       monitor_;

        std::vector<ScalarT> ws_;
        std::vector<IdxT>    ws_indices_;
      };
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit
