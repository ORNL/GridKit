/**
 * @file   CoSim.hpp
 * @author Philip Fackler (facklerpw@ornl.gov)
 *
 * @brief Declaration of a CoSim connector interface
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>

// Forward declarations
namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Connector
    {
      template <typename RealT, typename IdxT>
      struct CoSimData;
    }

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
    namespace Connector
    {
      /**
       * @brief Internal variables of a `CoSim`
       *
       * @note Technically these are not owned by this component, but they must
       * be classified as internal so that they can be written to signal nodes
       * in `ComponentSignals`
       */
      enum class CoSimInternalVariables : size_t
      {
        VREAL,
        VIMAG,
        MAXIMUM,
      };

      /// External variables of a `CoSim`
      enum class CoSimExternalVariables : size_t
      {
        IREAL,
        IIMAG,
        MAXIMUM,
      };

      template <class ScalarT, typename IdxT>
      class CoSim : public Component<ScalarT, IdxT>
      {
        using Component<ScalarT, IdxT>::gridkit_component_id_;
        using Component<ScalarT, IdxT>::size_;

      public:
        using RealT           = typename Component<ScalarT, IdxT>::RealT;
        using model_data_type = CoSimData<RealT, IdxT>;
        using signal_type     = SignalNode<ScalarT, IdxT>;
        using bus_type        = BusBase<ScalarT, IdxT>;

        CoSim(bus_type* bus);
        ~CoSim();

        int setGridKitComponentID(IdxT) override final;
        int allocate() override final;
        int verify() const override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int evaluateResidual() override final;
        int evaluateJacobian() override final;

        /// Get the `ComponentSignals` from this `CoSim`
        auto getSignals()
            -> ComponentSignals<ScalarT,
                                IdxT,
                                CoSimInternalVariables,
                                CoSimExternalVariables>&
        {
          return signals_;
        }

      private:
        // Placeholders for variable indices (see note in allocate() method)
        IdxT vr_index_{0};
        IdxT vi_index_{0};

        // Signal pointers
        signal_type* ir_signal_;
        signal_type* ii_signal_;
        bus_type*    bus_;

        /// Component signal extension
        ComponentSignals<ScalarT, IdxT, CoSimInternalVariables, CoSimExternalVariables> signals_;
      };

    } // namespace Connector
  } // namespace PhasorDynamics
} // namespace GridKit
