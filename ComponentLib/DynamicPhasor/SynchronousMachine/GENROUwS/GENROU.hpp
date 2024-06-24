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
        const ScalarT& Pm() const
        {
            return param_[0];
        }

        const ScalarT& Ef() const
        {
            return param_[1];
        }

        ScalarT Pg();
        ScalarT Qg();
        ScalarT frequencyPenalty(ScalarT omega);
        ScalarT frequencyPenaltyDer(ScalarT omega);

    private:
        //
        // Private inlined accessor methods
        //

        const ScalarT dotDelta() const
        {
            return yp_[0];
        }

        const ScalarT dotOmega() const
        {
            return yp_[1];
        }

        const ScalarT dotEqp() const
        {
            return yp_[2];
        }

        const ScalarT dotPsidp() const
        {
            return yp_[3];
        }

        const ScalarT dotPsiqp() const
        {
            return yp_[4];
        }

        const ScalarT dotEdp() const
        {
            return yp_[5];
        }

        const ScalarT dotPsiqpp() const
        {
            return yp_[6];
        }

        const ScalarT dotPsidpp() const
        {
            return yp_[7];
        }

        const ScalarT dotPsipp() const
        {
            return yp_[8];
        }

        const ScalarT dotksat() const
        {
            return yp_[9];
        }

        const ScalarT dotVd() const
        {
            return yp_[10];
        }

        const ScalarT dotVq() const
        {
            return yp_[11];
        }

        const ScalarT dotTelec() const
        {
            return yp_[12];
        }

        const ScalarT dotId() const
        {
            return yp_[13];
        }

        const ScalarT dotIq() const
        {
            return yp_[14];
        }

        const ScalarT dotVr() const
        {
            return yp_[15];
        }

        const ScalarT dotVi() const
        {
            return yp_[16];
        }

        const ScalarT dotIr() const
        {
            return yp_[17];
        }

        const ScalarT dotIi() const
        {
            return yp_[18];
        }

        const ScalarT dotPmech() const
        {
            return yp_[19];
        }

        //
        // y_ values
        //

        const ScalarT delta() const
        {
            return y_[0];
        }

        const ScalarT omega() const
        {
            return y_[1];
        }

        const ScalarT Eqp() const
        {
            return y_[2];
        }

        const ScalarT Psidp() const
        {
            return y_[3];
        }

        const ScalarT Psiqp() const
        {
            return y_[4];
        }

        const ScalarT Edp() const
        {
            return y_[5];
        }

        const ScalarT Psiqpp() const
        {
            return y_[6];
        }

        const ScalarT Psidpp() const
        {
            return y_[7];
        }

        const ScalarT Psipp() const
        {
            return y_[8];
        }

        const ScalarT ksat() const
        {
            return y_[9];
        }

        const ScalarT Vd() const
        {
            return y_[10];
        }

        const ScalarT Vq() const
        {
            return y_[11];
        }

        const ScalarT Telec() const
        {
            return y_[12];
        }

        const ScalarT Id() const
        {
            return y_[13];
        }

        const ScalarT Iq() const
        {
            return y_[14];
        }

        const ScalarT Vr() const
        {
            return y_[15];
        }

        const ScalarT Vi() const
        {
            return y_[16];
        }

        const ScalarT Ir() const
        {
            return y_[17];
        }

        const ScalarT Ii() const
        {
            return y_[18];
        }

        const ScalarT Pmech() const
        {
            return y_[19];
        }

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


/* Paul Moon 6/6/2024 */

/*#include "ModelEvaluatorImpl.hpp"

namespace ModelLib
{
    template <class ScalarT, typename IdxT> class BaseBus;
}

namespace ModelLib
{
    /**
     * @brief Generator
     * 
     * @tparam ScalarT - Scalar type
     * @tparam IdxT    - Matrix and vector index type
     */
    /*template  <class ScalarT, typename IdxT>
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

        typedef typename ModelEvaluatorImpl<ScalarT, IdxT>::real_type real_type;
        typedef BaseBus<ScalarT, IdxT> bus_type;

    public:
        GENROU(BaseBus<ScalarT, IdxT>* bus, ScalarT P0 = 1.0, ScalarT Q0 = 0.0);
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
        const ScalarT& Pm() const
        {
            return param_[0];
        }

        const ScalarT& Ef() const
        {
            return param_[1];
        }

        ScalarT Pg();
        ScalarT Qg();
        ScalarT frequencyPenalty(ScalarT omega);
        ScalarT frequencyPenaltyDer(ScalarT omega);

    private:
        //
        // Private inlined accessor methods
        //

        const ScalarT dotDelta() const
        {
            return yp_[0];
        }

        const ScalarT dotOmega() const
        {
            return yp_[1];
        }

        const ScalarT dotEdp() const
        {
            return yp_[2];
        }

        const ScalarT dotEqp() const
        {
            return yp_[3];
        }

        const ScalarT delta() const
        {
            return y_[0];
        }

        const ScalarT omega() const
        {
            return y_[1];
        }

        const ScalarT Edp() const
        {
            return y_[2];
        }

        const ScalarT Eqp() const
        {
            return y_[3];
        }

        const ScalarT Id() const
        {
            return y_[4];
        }

        const ScalarT Iq() const
        {
            return y_[5];
        }

    private:
        inline ScalarT frequencyPenalty(ScalarT omega);
        inline ScalarT frequencyPenaltyDer(ScalarT omega);

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

        real_type delta_; ///< Machine internal angle
        real_type omega_; ///<
        real_type Eqp_;   ///<
        real_type psidp_; ///<
        real_type psiqp_; ///<
        real_type Edp_;   ///<
        real_type psi_;   ///<
        real_type ksat_;  ///< Saturation coefficient
        real_type Vd_;    ///<
        real_type Vq_;    ///<
        real_type Telec_; ///<
        real_type Ir_;    ///<
        real_type Ii_;    ///<
        real_type Vr_;    ///<
        real_type Vi_;    ///<
        real_type Id_;    ///<
        real_type Iq_;    ///<
        real_type Pmech_; ///<
        real_type Efd_;   ///<




        real_type Ef_;    //
        real_type Pm_;
        real_type omega_s_;
        real_type omega_b_;
        real_type omega_up_;
        real_type omega_lo_;
        real_type c_;
        real_type beta_;

        ScalarT P0_;
        ScalarT Q0_;

        bus_type* bus_;
    };
}*/