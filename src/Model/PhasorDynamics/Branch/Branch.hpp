/**
 * @file Branch.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Declaration of a phasor dynamics branch model.
 *
 * The model uses Cartesian coordinates.
 *
 */
#pragma once

#include <Model/PhasorDynamics/Component.hpp>

// Forward declarations.
namespace GridKit
{
namespace PowerSystemData
{
    template <typename RealT, typename IdxT>
    struct BranchData;
}
} // namespace GridKit

namespace GridKit
{
namespace PhasorDynamics
{
    template <class ScalarT, typename IdxT>
    class Bus;
}
} // namespace GridKit

namespace GridKit
{
namespace PhasorDynamics
{
    /**
     * @brief Implementation of a pi-model branch between two buses.
     *
     * The model is implemented in Cartesian coordinates. Positive current
     * direction is into the busses.
     *
     */
    template <class ScalarT, typename IdxT>
    class Branch : public Component<ScalarT, IdxT>
    {
        using Component<ScalarT, IdxT>::size_;
        using Component<ScalarT, IdxT>::nnz_;
        using Component<ScalarT, IdxT>::time_;
        using Component<ScalarT, IdxT>::alpha_;
        using Component<ScalarT, IdxT>::y_;
        using Component<ScalarT, IdxT>::yp_;
        using Component<ScalarT, IdxT>::tag_;
        using Component<ScalarT, IdxT>::f_;
        using Component<ScalarT, IdxT>::g_;
        using Component<ScalarT, IdxT>::yB_;
        using Component<ScalarT, IdxT>::ypB_;
        using Component<ScalarT, IdxT>::fB_;
        using Component<ScalarT, IdxT>::gB_;
        using Component<ScalarT, IdxT>::param_;

        using bus_type   = BusBase<ScalarT, IdxT>;
        using real_type  = typename Component<ScalarT, IdxT>::real_type;
        using BranchData = GridKit::PowerSystemData::BranchData<real_type, IdxT>;

    public:
        Branch(bus_type* bus1, bus_type* bus2);
        Branch(bus_type* bus1, bus_type* bus2, real_type R, real_type X, real_type G, real_type B);
        Branch(bus_type* bus1, bus_type* bus2, BranchData& data);
        virtual ~Branch();

        virtual int allocate() override;
        virtual int initialize() override;
        virtual int tagDifferentiable() override;
        virtual int evaluateResidual() override;
        virtual int evaluateJacobian() override;
        virtual int evaluateIntegrand() override;

        virtual int initializeAdjoint() override;
        virtual int evaluateAdjointResidual() override;
        // virtual int evaluateAdjointJacobian() override;
        virtual int evaluateAdjointIntegrand() override;

        virtual void updateTime(real_type /* t */, real_type /* a */) override
        {
        }

    public:
        void setR(real_type R)
        {
            R_ = R;
        }

        void setX(real_type X)
        {
            // std::cout << "Setting X ...\n";
            X_ = X;
        }

        void setG(real_type G)
        {
            G_ = G;
        }

        void setB(real_type B)
        {
            B_ = B;
        }

    private:
        ScalarT& Vr1()
        {
            return bus1_->Vr();
        }

        ScalarT& Vi1()
        {
            return bus1_->Vi();
        }

        ScalarT& Ir1()
        {
            return bus1_->Ir();
        }

        ScalarT& Ii1()
        {
            return bus1_->Ii();
        }

        ScalarT& Vr2()
        {
            return bus2_->Vr();
        }

        ScalarT& Vi2()
        {
            return bus2_->Vi();
        }

        ScalarT& Ir2()
        {
            return bus2_->Ir();
        }

        ScalarT& Ii2()
        {
            return bus2_->Ii();
        }

    private:
        bus_type*  bus1_;
        bus_type*  bus2_;
        real_type  R_;
        real_type  X_;
        real_type  G_;
        real_type  B_;
        const IdxT bus1ID_;
        const IdxT bus2ID_;
    };

} // namespace PhasorDynamics
} // namespace GridKit
