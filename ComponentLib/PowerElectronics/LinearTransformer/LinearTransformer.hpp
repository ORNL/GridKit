

#ifndef _LTRANS_HPP_
#define _LTRANS_HPP_

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
     * @brief Declaration of a LinearTransformer class.
     *
     */
    template  <class ScalarT, typename IdxT>
    class LinearTransformer : public CircuitComponent<ScalarT, IdxT>
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
        LinearTransformer(IdxT id, ScalarT L1, ScalarT L2, ScalarT R1, ScalarT R2, ScalarT M);
        virtual ~LinearTransformer();

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
        ScalarT L1_;
        ScalarT L2_;
        ScalarT R1_;
        ScalarT R2_;
        ScalarT M_;
    };
}

#endif
