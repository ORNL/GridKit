

#ifndef _IMOTOR_HPP_
#define _IMOTOR_HPP_

#include <ModelEvaluatorImpl.hpp>
#include <PowerSystemData.hpp>
#include <ComponentLib/PowerElectronicsComponents/CircuitComponent.hpp>


namespace ModelLib
{
    template <class ScalarT, typename IdxT> class BaseBus;
}


namespace ModelLib
{
    /*!
     * @brief Declaration of a passive InductionMotor class.
     *
     */
    template  <class ScalarT, typename IdxT>
    class InductionMotor : public CircuitComponent<ScalarT, IdxT>
    {
        using CircuitComponent<ScalarT, IdxT>::size_;
        using CircuitComponent<ScalarT, IdxT>::nnz_;
        using CircuitComponent<ScalarT, IdxT>::time_;
        using CircuitComponent<ScalarT, IdxT>::alpha_;
        using CircuitComponent<ScalarT, IdxT>::y_;
        using CircuitComponent<ScalarT, IdxT>::yp_;
        using CircuitComponent<ScalarT, IdxT>::tag_;
        using CircuitComponent<ScalarT, IdxT>::f_;
        using CircuitComponent<ScalarT, IdxT>::g_;
        using CircuitComponent<ScalarT, IdxT>::yB_;
        using CircuitComponent<ScalarT, IdxT>::ypB_;
        using CircuitComponent<ScalarT, IdxT>::fB_;
        using CircuitComponent<ScalarT, IdxT>::gB_;
        using CircuitComponent<ScalarT, IdxT>::param_;
        using CircuitComponent<ScalarT, IdxT>::idc_;

    public:
        InductionMotor(IdxT id, ScalarT Lls, ScalarT Rs, ScalarT Llr, ScalarT Rr, ScalarT Lms, ScalarT J, ScalarT P);
        virtual ~InductionMotor();

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
        
    private:
		ScalarT Lls_;
		ScalarT Rs_;
		ScalarT Llr_;
		ScalarT Rr_;
		ScalarT Lms_;
		ScalarT J_;
		ScalarT P_;
    };
}

#endif
