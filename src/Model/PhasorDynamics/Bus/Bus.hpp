
#pragma once

#include <Model/PhasorDynamics/BusBase.hpp>

// Forward declaration of BusData structure
namespace GridKit
{
namespace PowerSystemData
{
    template <typename RealT, typename IdxT>
    struct BusData;
}
} // namespace GridKit

namespace GridKit
{
namespace PhasorDynamics
{
    /*!
     * @brief Implementation of a PQ bus.
     *
     * Voltage _V_ and phase _theta_ are variables in PQ bus model.
     * Active and reactive power, _P_ and _Q_, are residual components.
     *
     *
     */
    template <class ScalarT, typename IdxT>
    class Bus : public BusBase<ScalarT, IdxT>
    {
        using BusBase<ScalarT, IdxT>::size_;
        using BusBase<ScalarT, IdxT>::y_;
        using BusBase<ScalarT, IdxT>::yp_;
        using BusBase<ScalarT, IdxT>::yB_;
        using BusBase<ScalarT, IdxT>::ypB_;
        using BusBase<ScalarT, IdxT>::f_;
        using BusBase<ScalarT, IdxT>::fB_;
        using BusBase<ScalarT, IdxT>::tag_;

    public:
        using real_type = typename BusBase<ScalarT, IdxT>::real_type;
        using BusData   = GridKit::PowerSystemData::BusData<real_type, IdxT>;

        Bus();
        Bus(ScalarT Vr, ScalarT Vi);
        Bus(BusData& data);
        virtual ~Bus();

        virtual int allocate() override;
        virtual int tagDifferentiable() override;
        virtual int initialize() override;
        virtual int evaluateResidual() override;
        virtual int evaluateIntegrand() override;
        virtual int evaluateJacobian() override;
        virtual int initializeAdjoint() override;
        virtual int evaluateAdjointIntegrand() override;
        virtual int evaluateAdjointResidual() override;

        virtual int BusType() const override
        {
            return BusBase<ScalarT, IdxT>::BusType::DEFAULT;
        }

        virtual ScalarT& Vr() override
        {
            return y_[0];
        }

        virtual const ScalarT& Vr() const override
        {
            return y_[0];
        }

        virtual ScalarT& Vi() override
        {
            return y_[1];
        }

        virtual const ScalarT& Vi() const override
        {
            return y_[1];
        }

        virtual ScalarT& Ir() override
        {
            return f_[0];
        }

        virtual const ScalarT& Ir() const override
        {
            return f_[0];
        }

        virtual ScalarT& Ii() override
        {
            return f_[1];
        }

        virtual const ScalarT& Ii() const override
        {
            return f_[1];
        }

        // virtual ScalarT& VrB() override
        // {
        //     return yB_[0];
        // }

        // virtual const ScalarT& VrB() const override
        // {
        //     return yB_[0];
        // }

        // virtual ScalarT& ViB() override
        // {
        //     return yB_[1];
        // }

        // virtual const ScalarT& ViB() const override
        // {
        //     return yB_[1];
        // }

        // virtual ScalarT& IrB() override
        // {
        //     return fB_[0];
        // }

        // virtual const ScalarT& IrB() const override
        // {
        //     return fB_[0];
        // }

        // virtual ScalarT& IiB() override
        // {
        //     return fB_[1];
        // }

        // virtual const ScalarT& IiB() const override
        // {
        //     return fB_[1];
        // }

    private:
        ScalarT Vr0_{0.0};
        ScalarT Vi0_{0.0};
    };

} // namespace PhasorDynamics
} // namespace GridKit
