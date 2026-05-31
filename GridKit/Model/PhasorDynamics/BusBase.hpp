#pragma once

#include <memory>
#include <set>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Constants.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusData.hpp>
#include <GridKit/Model/PhasorDynamics/GridElement.hpp>
#include <GridKit/Model/VariableMonitor.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using Log = ::GridKit::Utilities::Logger;

    /*!
     * @brief BusBase model implementation base class.
     *
     */
    template <typename ScalarT, typename IdxT>
    class BusBase : public GridElement<ScalarT, IdxT>
    {
    public:
      using RealT    = typename GridElement<ScalarT, IdxT>::RealT;
      using MatrixT  = typename GridElement<ScalarT, IdxT>::MatrixT;
      using BusTypeT = typename BusData<RealT, IdxT>::BusType;
      using MonitorT = Model::VariableMonitor<BusBase, BusData>;

      BusBase() = default;

      virtual ~BusBase() = default;

      int verify() const override
      {
        return 0;
      }

      /// Pure virtual function, returns bus type (DEFAULT or SLACK).
      virtual BusTypeT BusType() const
      {
        return BusTypeT::DEFAULT;
      }

      bool hasJacobian() override
      {
        return false;
      }

      void updateTime(RealT /* t */, RealT /* a */) override
      {
        // No time to update in bus models
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

      const Model::VariableMonitorBase* getMonitor() const override
      {
        return monitor_.get();
      }

    protected:
      IdxT bus_id_{INVALID_INDEX<IdxT>};

      /// Variable monitor
      std::unique_ptr<MonitorT> monitor_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
