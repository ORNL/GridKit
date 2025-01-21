

#ifndef _TRANLINE_HPP_
#define _TRANLINE_HPP_

#include <ModelEvaluatorImpl.hpp>
#include <PowerSystemData.hpp>
#include <ComponentLib/PowerElectronics/CircuitComponent.hpp>


namespace ModelLib
{
    template <class ScalarT, typename IdxT> class BaseBus;
}


namespace ModelLib
{
    /*!
     * @brief Declaration of a TransmissionLine class.
     *
     * Model from Adam Birchfield paper (medium distances < 2km).
     * See also textbooks "Power System Analysis" by Grainger and "Power System Dynamics and Stability" by Sauer & Pai
     * 
     * @note Not used in the Microgrid model.
     */
    template  <class ScalarT, typename IdxT>
    class TransmissionLine : public CircuitComponent<ScalarT, IdxT>
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
        using CircuitComponent<ScalarT, IdxT>::J_;
        using CircuitComponent<ScalarT, IdxT>::param_;
        using CircuitComponent<ScalarT, IdxT>::idc_;

        using CircuitComponent<ScalarT, IdxT>::extern_indices_;
        using CircuitComponent<ScalarT, IdxT>::n_extern_;
        using CircuitComponent<ScalarT, IdxT>::n_intern_;

    public:
        TransmissionLine(IdxT id, ScalarT R, ScalarT X, ScalarT B);
        virtual ~TransmissionLine();

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
        ScalarT X_;
        ScalarT B_;
        ScalarT YReMat_;
        ScalarT YImMatDi_;
        ScalarT YImMatOff_;
    };
}

#endif
