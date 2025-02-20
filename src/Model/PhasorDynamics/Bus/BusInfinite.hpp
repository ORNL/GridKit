
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
}

namespace GridKit
{
namespace PhasorDynamics
{
    /*!
     * @brief Implementation of an "infinite" bus.
     *
     *
     *
     */
    template  <class ScalarT, typename IdxT>
    class BusInfinite : public BusBase<ScalarT, IdxT>
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

        BusInfinite();
        BusInfinite(ScalarT Vr, ScalarT Vi);
        BusInfinite(BusData& data);
        virtual ~BusInfinite();

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
            return BusBase<ScalarT, IdxT>::BusType::SLACK;
        }

        virtual ScalarT& Vr() override
        {
            return Vr_;
        }

        virtual const ScalarT& Vr() const override
        {
            return Vr_;
        }

        virtual ScalarT& Vi() override
        {
            return Vi_;
        }

        virtual const ScalarT& Vi() const override
        {
            return Vi_;
        }

        virtual ScalarT& Ir() override
        {
            return Ir_;
        }

        virtual const ScalarT& Ir() const override
        {
            return Ir_;
        }

        virtual ScalarT& Ii() override
        {
            return Ii_;
        }

        virtual const ScalarT& Ii() const override
        {
            return Ii_;
        }

        // virtual ScalarT& VrB() override
        // {
        //     return VrB_;
        // }

        // virtual const ScalarT& VrB() const override
        // {
        //     return VrB_;
        // }

        // virtual ScalarT& ViB() override
        // {
        //     return ViB_;
        // }

        // virtual const ScalarT& ViB() const override
        // {
        //     return ViB_;
        // }

        // virtual ScalarT& IrB() override
        // {
        //     return IrB_;
        // }

        // virtual const ScalarT& IrB() const override
        // {
        //     return IrB_;
        // }

        // virtual ScalarT& IiB() override
        // {
        //     return IiB_;
        // }

        // virtual const ScalarT& IiB() const override
        // {
        //     return IiB_;
        // }

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

} // PhasorDynamics
} // namespace GridKit
