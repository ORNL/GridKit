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
    template <typename real_type, typename index_type>
    struct BusToSignalAdapterData;

    template <typename scalar_type, typename index_type>
    class BusBase;

    template <typename scalar_type, typename index_type>
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

    template <typename scalar_type, typename index_type>
    class BusToSignalAdapter : public Component<scalar_type, index_type>
    {
      using Component<scalar_type, index_type>::gridkit_component_id_;
      using Component<scalar_type, index_type>::size_;

    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename Component<ScalarT, IdxT>::RealT;
      using ModelDataT = BusToSignalAdapterData<RealT, IdxT>;
      using SignalT    = SignalNode<ScalarT, IdxT>;
      using BusT       = BusBase<ScalarT, IdxT>;

      BusToSignalAdapter(BusT* bus);
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
      SignalT* ir_signal_;
      SignalT* ii_signal_;
      BusT*    bus_;

      /// Component signal extension
      ComponentSignals<ScalarT, IdxT, BusToSignalAdapterInternalVariables, BusToSignalAdapterExternalVariables> signals_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
