// Made by Paul Moon 6/7/2024

#ifndef _GENROU_H_
#define _GENROU_H_

#include <ModelEvaluatorImpl.hpp>

namespace ModelLib
{
    template <class ScalarT, typename IdxT> class BaseBus;
}

namespace ModelLib
{
    /*!
     * @brief Implementation of a fourth order generator model.
     *
     */
    template  <class ScalarT, typename IdxT>
    class GENROU : public ModelEvaluatorImpl<ScalarT, IdxT>
    {
        using ModelEvaluatorImpl<ScalarT, IdxT>::size_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::nnz_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::time_;
        using ModelEvaluatorImpl<ScalarT, IdxT>::alpha_;
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
        using ModelEvaluatorImpl<ScalarT, IdxT>::J_;

        typedef typename ModelEvaluatorImpl<ScalarT, IdxT>::real_type real_type;
        typedef BaseBus<ScalarT, IdxT> bus_type;

    public:
        GENROU(BaseBus<ScalarT, IdxT>* bus, ScalarT P0 = 0.999, ScalarT Q0 = 0.5723);
        virtual ~GENROU();

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

        void updateTime(real_type t, real_type a)
        {
            time_ = t;
            alpha_ = a;
        }

        // Inline accesor functions
        ScalarT& V()
        {
            return bus_->V();
        }

        const ScalarT& V() const
        {
            return bus_->V();
        }

        ScalarT& theta()
        {
            return bus_->theta();
        }

        const ScalarT& theta() const
        {
            return bus_->theta();
        }

        ScalarT& P()
        {
            return bus_->P();
        }

        const ScalarT& P() const
        {
            return bus_->P();
        }

        ScalarT& Q()
        {
            return bus_->Q();
        }

        const ScalarT& Q() const
        {
            return bus_->Q();
        }


    private:
        ScalarT Pg();
        ScalarT Qg();
        ScalarT frequencyPenalty(ScalarT omega);
        ScalarT frequencyPenaltyDer(ScalarT omega);

    private:
        real_type omega0_;///< Nominal frequency (assume 2pi60 rad/s)
        real_type H_;     ///< Inertia constant [s]
        real_type D_;     ///< Damping constant [pu]
        real_type Rs_;    ///< stator resistance [pu]
        real_type Xl_;    ///< stator leakage reactance [pu]
        real_type Xd_;    ///< direct axis synchronous reactance [pu]
        real_type Xdp_;   ///< direct axis transient reactance [pu]
        real_type Xdpp_;  ///< direct axis sub-transient reactance [pu]
        real_type Xq_;    ///< quadrature axis synchronous reactance [pu]
        real_type Xqp_;   ///< quadrature axis transient reactance [pu]
        real_type Xqpp_;  ///< quadrature axis sub-transient reactance [pu]
        real_type Td0p_;  ///< open circuit direct axis transient time constant [s]
        real_type Td0pp_; ///< open circuit direct axis sub-transient time constant [s]
        real_type Tq0p_;  ///< open circuit quadrature axis transient time constant [s]
        real_type Tq0pp_; ///< open circuit quadrature axis sub-transient time constant [s]
        real_type S1_;    ///< saturation at 1.0 pu flux
        real_type S12_;   ///< saturation at 1.2 pu flux
        real_type X_;

        real_type Ef_;
        real_type Pm_;
        real_type omega_s_;
        real_type omega_up_;
        real_type omega_lo_;
        real_type c_;
        real_type beta_;

        ScalarT P0_;
        ScalarT Q0_;

        bus_type* bus_;
    };

} // namespace ModelLib


#endif // _GENROU_H_
