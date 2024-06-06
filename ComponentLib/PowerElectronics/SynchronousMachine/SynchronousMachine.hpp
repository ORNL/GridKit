

#ifndef _SYNMACH_HPP_
#define _SYNMACH_HPP_

#include <ModelEvaluatorImpl.hpp>
#include <PowerSystemData.hpp>
#include <ComponentLib/PowerElectronics/CircuitComponent.hpp>
#include <tuple>


namespace ModelLib
{
    template <class ScalarT, typename IdxT> class BaseBus;
}


namespace ModelLib
{
    /*!
     * @brief Declaration of a SynchronousMachine class.
     *
     */
    template  <class ScalarT, typename IdxT>
    class SynchronousMachine : public CircuitComponent<ScalarT, IdxT>
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
        SynchronousMachine(IdxT id, ScalarT Lls, std::tuple<ScalarT,ScalarT> Llkq, ScalarT Llfd, ScalarT Llkd, ScalarT Lmq, ScalarT Lmd, ScalarT Rs, std::tuple<ScalarT,ScalarT> Rkq, ScalarT Rfd, ScalarT Rkd, ScalarT RJ, ScalarT P, ScalarT mub);
        virtual ~SynchronousMachine();

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
        std::tuple<ScalarT,ScalarT> Llkq_;
        ScalarT Llfd_;
        ScalarT Llkd_;
        ScalarT Lmq_;
        ScalarT Lmd_;
        ScalarT Rs_;
        std::tuple<ScalarT,ScalarT> Rkq_;
        ScalarT Rfd_;
        ScalarT Rkd_;
        ScalarT RJ_;
        ScalarT P_;
        ScalarT mub_;
    };
}

#endif
