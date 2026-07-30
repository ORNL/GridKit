#pragma once

#include <memory>

#include <GridKit/Constants.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusData.hpp>
#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Internal variables for buses
    enum class BusInternalVariables : size_t
    {
      VR,
      VI,
      MAXIMUM,
    };

    /// No external variables for buses
    enum class BusExternalVariables : size_t
    {
      MAXIMUM,
    };

    /*!
     * @brief BusBase model implementation base class.
     *
     * Buses are components that own the network voltage variables and the
     * current balance residual rows. This class adds the bus-specific
     * interface on top of Component; storage, binding, and index-map
     * machinery are inherited.
     */
    template <typename scalar_type, typename index_type>
    class BusBase : public Component<scalar_type, index_type>
    {
    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename Component<ScalarT, IdxT>::RealT;
      using CsrMatrixT = typename Component<ScalarT, IdxT>::CsrMatrixT;
      using CooMatrixT = typename Component<ScalarT, IdxT>::CooMatrixT;
      using VectorT    = typename Component<ScalarT, IdxT>::VectorT;
      using BusTypeT   = typename BusData<RealT, IdxT>::BusType;
      using MonitorT   = Model::VariableMonitor<BusBase, BusData>;

      BusBase() = default;

      virtual ~BusBase();

      virtual int verify() const override
      {
        return 0;
      }

      int setGridKitComponentID(IdxT component_id) override
      {
        this->gridkit_component_id_ = component_id;
        return 0;
      }

      /// Pure virtual function, returns bus type (DEFAULT or SLACK).
      virtual BusTypeT BusType() const
      {
        return BusTypeT::DEFAULT;
      }

      virtual ScalarT&       Vr()       = 0;
      virtual const ScalarT& Vr() const = 0;
      virtual ScalarT&       Vi()       = 0;
      virtual const ScalarT& Vi() const = 0;
      virtual ScalarT&       Ir()       = 0;
      virtual const ScalarT& Ir() const = 0;
      virtual ScalarT&       Ii()       = 0;
      virtual const ScalarT& Ii() const = 0;

      virtual int setBusID(IdxT) = 0;

      virtual const IdxT busID() const
      {
        return bus_id_;
      }

      const Model::VariableMonitorBase* getMonitor() const override;

      /// Get the `ComponentSignals` from this bus
      auto getSignals()
          -> ComponentSignals<ScalarT, IdxT, BusInternalVariables, BusExternalVariables>&
      {
        return signals_;
      }

    protected:
      IdxT bus_id_{INVALID_INDEX<IdxT>};

      /// Variable monitor
      std::unique_ptr<MonitorT> monitor_;

      /// Component signals
      ComponentSignals<ScalarT, IdxT, BusInternalVariables, BusExternalVariables> signals_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
