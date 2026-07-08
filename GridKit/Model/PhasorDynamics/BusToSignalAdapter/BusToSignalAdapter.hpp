/**
 * @file   BusToSignalAdapter.hpp
 * @author Philip Fackler (facklerpw@ornl.gov)
 *
 * @brief Declaration of a BusToSignalAdapter component
 */

#pragma once

#include <GridKit/Constants.hpp>
#include <GridKit/Model/PhasorDynamics/BusToSignalAdapter/BusToSignalAdapterData.hpp>
#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/IOPorts.hpp>

// Forward declarations
namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename scalar_type, typename index_type>
    class BusBase;

    template <typename scalar_type, typename index_type>
    class SignalNode;

    template <typename signal_node_type>
    class SignalNodeSet;
  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Internal variables of a `BusToSignalAdapter`
    enum class BusToSignalAdapterInternalVariables : size_t
    {
      MAXIMUM,
    };

    /// External variables of a `BusToSignalAdapter`
    enum class BusToSignalAdapterExternalVariables : size_t
    {
      MAXIMUM,
    };

    template <typename scalar_type, typename index_type>
    class BusToSignalAdapter : public Component<scalar_type, index_type>
    {
      using Component<scalar_type, index_type>::gridkit_component_id_;
      using Component<scalar_type, index_type>::size_;

    public:
      using ScalarT        = scalar_type;
      using IdxT           = index_type;
      using RealT          = typename Component<ScalarT, IdxT>::RealT;
      using ModelDataT     = BusToSignalAdapterData<RealT, IdxT>;
      using BusT           = BusBase<ScalarT, IdxT>;
      using SignalNodeT    = SignalNode<ScalarT, IdxT>;
      using SignalNodeSetT = SignalNodeSet<SignalNodeT>;
      using IOPortsT       = IOPorts<ScalarT, ModelDataT>;

      BusToSignalAdapter(BusT* bus);
      BusToSignalAdapter(BusT*             bus,
                         const ModelDataT& data,
                         SignalNodeSetT&   signal_nodes);
      ~BusToSignalAdapter();

      int setGridKitComponentID(IdxT) override final;
      int allocate() override final;
      int verify() const override final;
      int initialize() override final;
      int tagDifferentiable() override final;
      int setAbsoluteTolerance(RealT rel_tol) override final;
      int evaluateResidual() override final;
      int evaluateJacobian() override final;

      IOPortsT& getPorts()
      {
        return ports_;
      }

    private:
      // Placeholders for variable indices (see note in allocate() method)
      IdxT vr_index_{INVALID_INDEX<IdxT>};
      IdxT vi_index_{INVALID_INDEX<IdxT>};

      // Bus pointer
      BusT* bus_;

      // Component ports
      IOPortsT ports_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
