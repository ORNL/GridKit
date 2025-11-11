

#pragma once

#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>

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
  };
} // namespace GridKit

namespace GridKit
{
  /*!
   * @brief Declaration of a DistributedGenerator class.
   *
   */
  template <class ScalarT, typename IdxT>
  class DistributedGenerator : public CircuitComponent<ScalarT, IdxT>
  {
    using RealT = typename CircuitComponent<ScalarT, IdxT>::RealT;

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

  public:
    DistributedGenerator(IdxT                                        id,
                         DistributedGeneratorParameters<RealT, IdxT> parm,
                         bool                                        reference_frame);
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
