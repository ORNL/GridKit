
#pragma once

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/ConnectedElement.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename, typename>
    class BusInfinite;

    template <typename ScalarP, typename IdxP>
    struct ConnectedElementTraits<BusInfinite<ScalarP, IdxP>>
    {
      using BusT = BusInfinite<ScalarP, IdxP>;

      using ElementT           = BusT;
      using ScalarT            = ScalarP;
      using IdxT               = IdxP;
      using RealT              = typename ScalarTraits<ScalarT>::RealT;
      using ModelDataT         = BusData<RealT, IdxT>;
      using InternalVariablesT = BusInternalVariables;
      using ExternalVariablesT = BusExternalVariables;
      using InterfaceT         = BusBase<ScalarT, IdxT>;
    };

    /*!
     * @brief Implementation of an "infinite" bus.
     *
     *
     *
     */
    template <typename ScalarP, typename IdxP>
    class BusInfinite : public ConnectedElement<BusInfinite<ScalarP, IdxP>>
    {
      using ConnectedElement<BusInfinite>::bus_id_;
      using ConnectedElement<BusInfinite>::size_;
      using ConnectedElement<BusInfinite>::y_;
      using ConnectedElement<BusInfinite>::yp_;
      using ConnectedElement<BusInfinite>::f_;
      using ConnectedElement<BusInfinite>::J_;
      using ConnectedElement<BusInfinite>::variable_indices_;
      using ConnectedElement<BusInfinite>::residual_indices_;
      using ConnectedElement<BusInfinite>::monitor_;

    public:
      using ScalarT    = typename ConnectedElement<BusInfinite>::ScalarT;
      using IdxT       = typename ConnectedElement<BusInfinite>::IdxT;
      using RealT      = typename ConnectedElement<BusInfinite>::RealT;
      using ModelDataT = BusData<RealT, IdxT>;
      using BusTypeT   = typename ConnectedElement<BusInfinite>::BusTypeT;

      BusInfinite();
      BusInfinite(ScalarT Vr, ScalarT Vi);
      BusInfinite(const ModelDataT& data);
      virtual ~BusInfinite();

      virtual int setBusID(IdxT) override final;
      virtual int allocate() override final;
      virtual int tagDifferentiable() override final;
      virtual int initialize() override final;
      virtual int evaluateResidual() override final;
      virtual int evaluateJacobian() override final;

      virtual BusTypeT BusType() const override final
      {
        return BusTypeT::SLACK;
      }

      virtual ScalarT& Vr() override final
      {
        return Vr_;
      }

      virtual const ScalarT& Vr() const override final
      {
        return Vr_;
      }

      virtual ScalarT& Vi() override final
      {
        return Vi_;
      }

      virtual const ScalarT& Vi() const override final
      {
        return Vi_;
      }

      virtual ScalarT& Ir() override final
      {
        return Ir_;
      }

      virtual const ScalarT& Ir() const override final
      {
        return Ir_;
      }

      virtual ScalarT& Ii() override final
      {
        return Ii_;
      }

      virtual const ScalarT& Ii() const override final
      {
        return Ii_;
      }

    private:
      ScalarT Vr_{0.0};
      ScalarT Vi_{0.0};
      ScalarT Ir_{0.0};
      ScalarT Ii_{0.0};

      ScalarT VrB_{0.0};
      ScalarT ViB_{0.0};
      ScalarT IrB_{0.0};
      ScalarT IiB_{0.0};
    };

  } // namespace PhasorDynamics
} // namespace GridKit
