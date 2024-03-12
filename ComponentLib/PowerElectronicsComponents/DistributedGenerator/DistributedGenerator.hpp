

#ifndef _CAP_HPP_
#define _CAP_HPP_

#include <ModelEvaluatorImpl.hpp>
#include <PowerSystemData.hpp>
#include <ComponentLib/PowerElectronicsComponents/CircuitComponent.hpp>


namespace ModelLib
{
    template <class ScalarT, typename IdxT> class BaseBus;

    template <class ScalarT, typename IdxT>
    struct DistributedGeneratorParameters
    {
        ScalarT wb;
        ScalarT wc;
        ScalarT mp;
        ScalarT Vn;
        ScalarT nq;
        ScalarT F;
        ScalarT Kiv;
        ScalarT Kpv;
        ScalarT Kic;
        ScalarT Kpc;
        ScalarT Cf;
        ScalarT rLf;
        ScalarT Lf;
        ScalarT rLc;
        ScalarT Lc;
    };
}


namespace ModelLib
{
    /*!
     * @brief Declaration of a passive DistributedGenerator class.
     *
     */
    template  <class ScalarT, typename IdxT>
    class DistributedGenerator : public CircuitComponent<ScalarT, IdxT>
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


    public:
        DistributedGenerator(IdxT id, DistributedGeneratorParameters<ScalarT,IdxT> parm, bool reference_frame);
        virtual ~DistributedGenerator();

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
        ScalarT wb_;
        ScalarT wc_;
        ScalarT mp_;
        ScalarT Vn_;
        ScalarT nq_;
        ScalarT F_;
        ScalarT Kiv_;
        ScalarT Kpv_;
        ScalarT Kic_;
        ScalarT Kpc_;
        ScalarT Cf_;
        ScalarT rLf_;
        ScalarT Lf_;
        ScalarT rLc_;
        ScalarT Lc_;
        bool refframe_;
    };
}

#endif
