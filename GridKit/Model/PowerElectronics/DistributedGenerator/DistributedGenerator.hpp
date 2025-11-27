

#pragma once

#include <GridKit/Model/EvaluatorMixins.hpp>
#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>
#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

namespace GridKit
{
  template <class ScalarT, typename IdxT>
  class BaseBus;

  template <typename RealT, typename IdxT>
  struct DistributedGeneratorParameters
  {
    RealT wb_;
    RealT wc_;
    RealT mp_;
    RealT Vn_;
    RealT nq_;
    RealT F_;
    RealT Kiv_;
    RealT Kpv_;
    RealT Kic_;
    RealT Kpc_;
    RealT Cf_;
    RealT rLf_;
    RealT Lf_;
    RealT rLc_;
    RealT Lc_;

    DistributedGeneratorParameters()
    {
      wb_ = 0.0;
      wc_ = 0.0;
      mp_ = 0.0;
      Vn_ = 0.0;
      nq_ = 0.0;
      F_ = 0.0;
      Kiv_ = 0.0;
      Kpv_ = 0.0;
      Kic_ = 0.0;
      Kpc_ = 0.0;
      Cf_ = 0.0;
      rLf_ = 0.0;
      Lf_ = 0.0;
      rLc_ = 0.0;
      Lc_ = 0.0;
    }

    template <class S, typename I>
    DistributedGeneratorParameters(const DistributedGeneratorParameters<S, I>& other)
    {
      wb_ = other.wb_;
      wc_ = other.wc_;
      mp_ = other.mp_;
      Vn_ = other.Vn_;
      nq_ = other.nq_;
      F_ = other.F_;
      Kiv_ = other.Kiv_;
      Kpv_ = other.Kpv_;
      Kic_ = other.Kic_;
      Kpc_ = other.Kpc_;
      Cf_ = other.Cf_;
      rLf_ = other.rLf_;
      Lf_ = other.Lf_;
      rLc_ = other.rLc_;
      Lc_ =  other.Lc_;
    }
  };
} // namespace GridKit

namespace GridKit
{
  /*!
   * @brief Declaration of a DistributedGenerator class.
   *
   */
  template <class ScalarT, typename IdxT>
  class DistributedGenerator : public CircuitComponent<ScalarT, IdxT>, public Mixin::Evaluator::CsrJacobian<ScalarT, IdxT, DistributedGenerator>
  {
    using RealT   = typename CircuitComponent<ScalarT, IdxT>::RealT;
    using MatrixT = typename CircuitComponent<RealT, IdxT>::MatrixT;

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
    using CircuitComponent<ScalarT, IdxT>::jac_;
    using CircuitComponent<ScalarT, IdxT>::param_;
    using CircuitComponent<ScalarT, IdxT>::idc_;

    using CircuitComponent<ScalarT, IdxT>::extern_indices_;
    using CircuitComponent<ScalarT, IdxT>::n_extern_;
    using CircuitComponent<ScalarT, IdxT>::n_intern_;

    using typename Model::Evaluator<ScalarT, IdxT>::CsrJacobian;

  public:
    DistributedGenerator(IdxT                                        id,
                         DistributedGeneratorParameters<RealT, IdxT> parm,
                         bool                                        reference_frame);

    DistributedGenerator(const DistributedGenerator<ScalarT, IdxT>& other);

    CircuitComponent<ScalarT, IdxT>* clone() const override
    {
      return new DistributedGenerator<ScalarT, IdxT>(*this);
    }
    
    CircuitComponent<ScalarT, IdxT>& operator=(const CircuitComponent<ScalarT, IdxT>& other);
    
    virtual ~DistributedGenerator();

    int allocate();
    int initialize();
    int tagDifferentiable();
    int evaluateResidual();
    int evaluateJacobian();
    int evaluateIntegrand();
    int initializeAdjoint();
    int evaluateAdjointResidual();
    // int evaluateAdjointJacobian();
    int evaluateAdjointIntegrand();

    /**
     * @brief copy out system parameters
     * 
     * @tparam ScalarT 
     * @tparam IdxT 
     * @return DistributedGeneratorParameters<ScalarT, IdxT> 
     */
    DistributedGeneratorParameters<ScalarT, IdxT>  getParameters() const
    {
      DistributedGeneratorParameters<ScalarT, IdxT> parms;
      parms.wb_  = wb_;
      parms.wc_  = wc_;
      parms.mp_  = mp_;
      parms.Vn_  = Vn_;
      parms.nq_  = nq_;
      parms.F_   = F_;
      parms.Kiv_ = Kiv_;
      parms.Kpv_ = Kpv_;
      parms.Kic_ = Kic_;
      parms.Kpc_ = Kpc_;
      parms.Cf_  = Cf_;
      parms.rLf_ = rLf_;
      parms.Lf_  = Lf_;
      parms.rLc_ = rLc_;
      parms.Lc_  = Lc_;

      return parms;
    }
  
  
    /**
     * @brief Set the Parameters object.
     * 
     * extra templates allow typcasting
     * 
     * @tparam S 
     * @tparam I 
     * @param params 
     */
    template <class S, typename I>
    void setParameters(const DistributedGeneratorParameters<S, I> params)
    {
      wb_ = params.wb_;
      wc_ = params.wc_;
      mp_ = params.mp_;
      Vn_ = params.Vn_;
      nq_ = params.nq_;
      F_ = params.F_;
      Kiv_ = params.Kiv_;
      Kpv_ = params.Kpv_;
      Kic_ = params.Kic_;
      Kpc_ = params.Kpc_;
      Cf_ = params.Cf_;
      rLf_ = params.rLf_;
      Lf_ = params.Lf_;
      rLc_ = params.rLc_;
      Lc_ =  params.Lc_;
    }

    bool isReferenceFrame() const;
    void setReferenceFrame(const bool ref);

    template <bool INCLUDE_DIAGONALS, bool KEEP_SORTED, bool USE_TEMPLATE>
    CsrJacobian buildCsrJacobian(LinearAlgebra::CsrBuilder<RealT, IdxT, INCLUDE_DIAGONALS, KEEP_SORTED, USE_TEMPLATE> builder);

    static const std::size_t SIZE = 16;


  private:
    RealT wb_;
    RealT wc_;
    RealT mp_;
    RealT Vn_;
    RealT nq_;
    RealT F_;
    RealT Kiv_;
    RealT Kpv_;
    RealT Kic_;
    RealT Kpc_;
    RealT Cf_;
    RealT rLf_;
    RealT Lf_;
    RealT rLc_;
    RealT Lc_;
    bool  refframe_;
  };
} // namespace GridKit
