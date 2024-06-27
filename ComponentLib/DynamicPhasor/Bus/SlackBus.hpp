/* Created by Paul Moon 6/26/2024 */
#ifndef _SLACK_BUS_HPP_
#define _SLACK_BUS_HPP_

#include <PowerSystemData.hpp>
#include "PhasorBus.hpp"

namespace ModelLib
{
    /*!
     * @brief Implementation of a slack bus.
     *
     * Slack bus sets voltage _V_ and phase _theta_ as constants.
     * Active and reactive power, _P_ and _Q_, are component model outputs,
     * but are computed outside the SlackBus class.
     *
     *
     */
    template  <class ScalarT, typename IdxT>
    class SlackBus : public PhasorBus<ScalarT, IdxT>
    {
        using PhasorBus<ScalarT, IdxT>::size_;
        using PhasorBus<ScalarT, IdxT>::y_;
        using PhasorBus<ScalarT, IdxT>::yp_;
        using PhasorBus<ScalarT, IdxT>::f_;
        using PhasorBus<ScalarT, IdxT>::g_;
        using PhasorBus<ScalarT, IdxT>::atol_;
        using PhasorBus<ScalarT, IdxT>::rtol_;

    public:
        using real_type = typename ModelEvaluatorImpl<ScalarT, IdxT>::real_type;
        using BusData = GridKit::PowerSystemData::BusData<real_type, IdxT>;

        SlackBus();
        SlackBus(ScalarT Vr, ScalarT Vi);
        SlackBus(BusData& data);
        virtual ~SlackBus();
        virtual int evaluateResidual();
        virtual int evaluateAdjointResidual();

        /// @todo Should slack bus allow changing voltage?
        virtual ScalarT& Vr()
        {
            return Vr_;
        }

        virtual const ScalarT& Vr() const
        {
            return Vr_;
        }

        /// @todo Should slack bus allow changing phase?
        virtual ScalarT& Vi()
        {
            return Vi_;
        }

        virtual const ScalarT& Vi() const
        {
            return Vi_;
        }

        virtual ScalarT& P()
        {
            return P_;
        }

        virtual const ScalarT& P() const
        {
            return P_;
        }

        virtual ScalarT& Q()
        {
            return Q_;
        }

        virtual const ScalarT& Q() const
        {
            return Q_;
        }

        /// @todo Should slack bus allow changing voltage?
        virtual ScalarT& lambdaP()
        {
            return thetaB_;
        }

        virtual const ScalarT& lambdaP() const
        {
            return thetaB_;
        }

        /// @todo Should slack bus allow changing phase?
        virtual ScalarT& lambdaQ()
        {
            return VB_;
        }

        virtual const ScalarT& lambdaQ() const
        {
            return VB_;
        }

        virtual ScalarT& PB()
        {
            return PB_;
        }

        virtual const ScalarT& PB() const
        {
            return PB_;
        }

        virtual ScalarT& QB()
        {
            return QB_;
        }

        virtual const ScalarT& QB() const
        {
            return QB_;
        }

        virtual const int BusType() const
        {
            return PhasorBus<ScalarT, IdxT>::BusType::Slack;
        }

    private:
        ScalarT Vr_;
        ScalarT Vi_;
        ScalarT P_;
        ScalarT Q_;

        ScalarT VB_;
        ScalarT thetaB_;
        ScalarT PB_;
        ScalarT QB_;

    }; // class SlackBus

} // namespace ModelLib


#endif // _SLACK_BUS_HPP_