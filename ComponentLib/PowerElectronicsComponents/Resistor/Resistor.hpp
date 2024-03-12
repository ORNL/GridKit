

#ifndef _RES_HPP_
#define _RES_HPP_

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
     * @brief Declaration of a passive Resistor class.
     *
     */
    template  <class ScalarT, typename IdxT>
    class Resistor : public CircuitComponent<ScalarT, IdxT>
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
        Resistor(IdxT id, ScalarT R);
        virtual ~Resistor();

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
		ScalarT R_;
    };
}

#endif
