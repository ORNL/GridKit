/**
 * @file   BusToSignalAdapter.hpp
 * @author Philip Fackler (facklerpw@ornl.gov)
 *
 * @brief Declaration of a BusToSignalAdapter component
 */

#pragma once

#include <GridKit/Constants.hpp>
#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>

// Forward declarations
namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename RealT, typename IdxT>
    struct BusToSignalAdapterData;

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
    /**
     * @brief Internal variables of a `BusToSignalAdapter`
     *
     * @note Technically these are not owned by this component, but they must
     * be classified as internal so that they can be written to signal nodes
     * in `ComponentSignals`
     */
    enum class BusToSignalAdapterInternalVariables : size_t
    {
      VREAL,
      VIMAG,
      MAXIMUM,
    };

    /// External variables of a `BusToSignalAdapter`
    enum class BusToSignalAdapterExternalVariables : size_t
    {
      IREAL,
      IIMAG,
      MAXIMUM,
    };

    template <class ScalarT, typename IdxT>
    class BusToSignalAdapter : public Component<ScalarT, IdxT>
    {
      using Component<ScalarT, IdxT>::gridkit_component_id_;
      using Component<ScalarT, IdxT>::size_;

    public:
      using RealT           = typename Component<ScalarT, IdxT>::RealT;
      using model_data_type = BusToSignalAdapterData<RealT, IdxT>;
      using signal_type     = SignalNode<ScalarT, IdxT>;
      using bus_type        = BusBase<ScalarT, IdxT>;

      BusToSignalAdapter(bus_type* bus);
      ~BusToSignalAdapter();

      int setGridKitComponentID(IdxT) override final;
      int allocate() override final;
      int verify() const override final;
      int initialize() override final;
      int tagDifferentiable() override final;
      int evaluateResidual() override final;
      int evaluateJacobian() override final;

      /// Get the `ComponentSignals` from this `BusToSignalAdapter`
      auto getSignals()
          -> ComponentSignals<ScalarT,
                              IdxT,
                              BusToSignalAdapterInternalVariables,
                              BusToSignalAdapterExternalVariables>&
      {
        return signals_;
      }

    private:
      // Placeholders for variable indices (see note in allocate() method)
      IdxT vr_index_{INVALID_INDEX<IdxT>};
      IdxT vi_index_{INVALID_INDEX<IdxT>};

      // Signal pointers
      signal_type* ir_signal_;
      signal_type* ii_signal_;
      bus_type*    bus_;

      /// Component signal extension
      ComponentSignals<ScalarT, IdxT, BusToSignalAdapterInternalVariables, BusToSignalAdapterExternalVariables> signals_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
