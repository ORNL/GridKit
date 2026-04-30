
#pragma once

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/ConnectedElement.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename, typename>
    class Bus;

    template <typename ScalarP, typename IdxP>
    struct ConnectedElementTraits<Bus<ScalarP, IdxP>>
    {
      using BusT = Bus<ScalarP, IdxP>;

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
     * @brief Implementation of a PQ bus.
     *
     * Voltage _V_ and phase _theta_ are variables in PQ bus model.
     * Active and reactive power, _P_ and _Q_, are residual components.
     *
     *
     */
    template <typename ScalarP, typename IdxP>
    class Bus : public ConnectedElement<Bus<ScalarP, IdxP>>
    {
      using ConnectedElement<Bus>::bus_id_;
      using ConnectedElement<Bus>::size_;
      using ConnectedElement<Bus>::y_;
      using ConnectedElement<Bus>::yp_;
      using ConnectedElement<Bus>::f_;
      using ConnectedElement<Bus>::J_;
      using ConnectedElement<Bus>::J_rows_buffer_;
      using ConnectedElement<Bus>::J_cols_buffer_;
      using ConnectedElement<Bus>::J_vals_buffer_;
      using ConnectedElement<Bus>::tag_;
      using ConnectedElement<Bus>::variable_indices_;
      using ConnectedElement<Bus>::residual_indices_;
      using ConnectedElement<Bus>::monitor_;

    public:
      using ScalarT    = typename ConnectedElement<Bus>::ScalarT;
      using IdxT       = typename ConnectedElement<Bus>::IdxT;
      using RealT      = typename ConnectedElement<Bus>::RealT;
      using ModelDataT = BusData<RealT, IdxT>;
      using BusTypeT   = typename ConnectedElement<Bus>::BusTypeT;

      Bus();
      Bus(ScalarT Vr, ScalarT Vi);
      Bus(const ModelDataT& data);
      virtual ~Bus();

      virtual int setBusID(IdxT) override final;
      virtual int allocate() override final;
      virtual int tagDifferentiable() override final;
      virtual int initialize() override final;
      virtual int evaluateResidual() override final;
      virtual int evaluateJacobian() override final;

      virtual BusTypeT BusType() const override final
      {
        return BusTypeT::DEFAULT;
      }

      virtual ScalarT& Vr() override final
      {
        return y_[0];
      }

      virtual const ScalarT& Vr() const override final
      {
        return y_[0];
      }

      virtual ScalarT& Vi() override final
      {
        return y_[1];
      }

      virtual const ScalarT& Vi() const override final
      {
        return y_[1];
      }

      virtual ScalarT& Ir() override final
      {
        return f_[0];
      }

      virtual const ScalarT& Ir() const override final
      {
        return f_[0];
      }

      virtual ScalarT& Ii() override final
      {
        return f_[1];
      }

      virtual const ScalarT& Ii() const override final
      {
        return f_[1];
      }

    private:
      ScalarT Vr0_{0.0};
      ScalarT Vi0_{0.0};
    };

  } // namespace PhasorDynamics
} // namespace GridKit
