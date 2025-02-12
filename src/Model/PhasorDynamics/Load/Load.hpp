

#pragma once

#include <Model/PhasorDynamics/Component.hpp>

// Forward declarations.
namespace GridKit
{
namespace PhasorDynamics
{
    template <class ScalarT, typename IdxT> class Bus;
}
}

namespace GridKit
{
namespace PhasorDynamics
{
    /*!
     * @brief Implementation of a constant load.
     *
     */
    template  <class ScalarT, typename IdxT>
    class Load : public Component<ScalarT, IdxT>
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
        using Component<ScalarT, IdxT>::component_id_;

        using bus_type   = BusBase<ScalarT, IdxT>;
        using real_type  = typename Component<ScalarT, IdxT>::real_type;

    public:
        Load(bus_type* bus);
        Load(bus_type* bus, real_type R, real_type X);
        Load(bus_type* bus, IdxT component_id);
        virtual ~Load();

        int allocate();
        int initialize();
        int tagDifferentiable();
        int evaluateResidual();
        int evaluateJacobian();
        int evaluateIntegrand();

        int initializeAdjoint();
        int evaluateAdjointResidual();
        //int evaluateAdjointJacobian();
        int evaluateAdjointIntegrand();

        void updateTime(real_type /* t */, real_type /* a */)
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

    private:
        ScalarT& Vr()
        {
            return bus_->Vr();
        }

        ScalarT& Vi()
        {
            return bus_->Vi();
        }

        ScalarT& Ir()
        {
            return bus_->Ir();
        }

        ScalarT& Ii()
        {
            return bus_->Ii();
        }


    private:
        bus_type* bus_{nullptr};
        real_type R_{0.1};
        real_type X_{0.01};
    };

} // namespace PhasorDynamics
} // namespace GridKit
