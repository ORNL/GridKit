/* Paul Moon 6/26/2024 */

#ifndef _PHASOR_BUS_HPP_
#define _PHASOR_BUS_HPP_

#include <ModelEvaluatorImpl.hpp>

namespace ModelLib
{
    /*!
     * @brief Base class for all dynamic phasor buses.
     *
     * Derived bus types:
     *   0 - swing bus (Vr and Vi are constants)
     *   1 - PV bus    (P and V are constants)
     *   2 - PQ bus    (P and Q are constants)
     *
     * @todo Consider static instead of dynamic polymorphism for
     * bus types. Create Bus class that takes template parameter
     * BusType.
     */
    template  <class ScalarT, typename IdxT>
    class PhasorBus : public ModelEvaluatorImpl<ScalarT, IdxT>
    {
    protected:
        using ModelEvaluatorImpl<ScalarT, IdxT>::size_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::nnz_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::time_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::alpha_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::rtol_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::atol_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::y_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::yp_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::tag_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::f_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::g_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::yB_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::ypB_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::fB_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::gB_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::param_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::param_up_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::param_lo_;

    public:
        typedef typename ModelEvaluatorImpl<ScalarT, IdxT>::real_type real_type;

        enum BusType{PQ=1, PV, Slack, Isolated};

        PhasorBus(IdxT id) : busID_(id) {}
        virtual ~PhasorBus(){}

        // Set defaults for ModelEvaluator methods
        virtual int allocate() { return 0;}
        virtual int initialize() { return 0;}
        virtual int tagDifferentiable() { return 0;}
        virtual int evaluateResidual() { return 0;}
        virtual int evaluateJacobian() { return 0;}
        virtual int evaluateIntegrand() { return 0;}

        virtual int initializeAdjoint() { return 0;}
        virtual int evaluateAdjointResidual() { return 0;}
        //virtual int evaluateAdjointJacobian() { return 0;}
        virtual int evaluateAdjointIntegrand() { return 0;}
        virtual void updateTime(real_type, real_type) {} // <- throw exception here

        // Pure virtual methods specific to Bus types
        virtual ScalarT& Vr() = 0;
        virtual const ScalarT& Vr() const = 0;
        virtual ScalarT& Vi() = 0;
        virtual const ScalarT& Vi() const = 0;
        virtual ScalarT& P() = 0;
        virtual const ScalarT& P() const = 0;
        virtual ScalarT& Q() = 0;
        virtual const ScalarT& Q() const = 0;

        virtual ScalarT& lambdaP() = 0;
        virtual const ScalarT& lambdaP() const = 0;
        virtual ScalarT& lambdaQ() = 0;
        virtual const ScalarT& lambdaQ() const = 0;
        virtual ScalarT& PB() = 0;
        virtual const ScalarT& PB() const = 0;
        virtual ScalarT& QB() = 0;
        virtual const ScalarT& QB() const = 0;

        virtual const int BusType() const = 0;

        virtual const IdxT BusID() const
        {
            return busID_;
        }

    protected:
        const IdxT busID_;
    }; // class PhasorBus

} // namespace ModelLib


#endif // _PHASOR_BUS_HPP_
