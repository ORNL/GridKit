#pragma once

#include <memory>
#include <set>
#include <vector>

#include <GridKit/Constants.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusData.hpp>
#include <GridKit/Model/PhasorDynamics/GridElement.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename, typename>
    class BusBase;

    enum class BusInternalVariables : size_t
    {
      MAXIMUM
    };

    enum class BusExternalVariables : size_t
    {
      MAXIMUM
    };

    /*!
     * @brief BusBase model implementation base class.
     *
     */
    template <typename ScalarP, typename IdxP>
    class BusBase : public GridElement<ScalarP, IdxP>
    {
    public:
      using ScalarT    = typename GridElement<ScalarP, IdxP>::ScalarT;
      using IdxT       = typename GridElement<ScalarP, IdxP>::IdxT;
      using RealT      = typename GridElement<ScalarP, IdxP>::RealT;
      using MatrixT    = typename GridElement<ScalarP, IdxP>::MatrixT;
      using ModelDataT = BusData<RealT, IdxT>;
      using BusTypeT   = typename BusData<RealT, IdxT>::BusType;

      BusBase() = default;

      explicit BusBase(const BusData<RealT, IdxT>& data);

      virtual ~BusBase();

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

    protected:
      IdxT bus_id_{INVALID_INDEX<IdxT>};
    };

  } // namespace PhasorDynamics
} // namespace GridKit
